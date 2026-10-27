# 수평 스캐닝 라이다 분석 코드 (PKNU 대기원격탐사연구실)

수평 방향으로 스캔한 에어로졸 라이다 관측으로부터 **소산계수(EC)**,
**Ångström exponent(AE)**, **편광비(Dp)** 를 산출하는 코드입니다.
**MATLAB · Python 두 가지 버전**이 같은 내용으로 들어 있습니다.

참조 소산계수를 여러 거리 구간에서 구해 이상치를 반복 제거하는
**multi-section 방법**을 사용합니다 — Shin et al. (2024), *Atmos. Meas. Tech.*

> **한시 공개 저장소**
> 제69회 한국대기환경학회 정기학술대회 교육프로그램
> 「라이다 시스템 개요와 데이터 해석 및 분석법 실습」(2026. 10. 28.) 실습용으로
> 공개한 저장소입니다. 교육 종료 후 비공개로 전환될 수 있습니다.
> 필요하신 분은 교육 기간 중 내려받아 두시기 바랍니다.

---

## 저장소 구성

| 파일 | 내용 |
|---|---|
| `Klett_Main_PKNU.m` | 메인 분석 코드 — MATLAB (단일 파일) |
| `Klett_Main_PKNU.py` | 메인 분석 코드 — Python (단일 파일, 동일 내용) |
| `testData.zip` | 예제 관측 자료 — 압축을 풀어 `testData` 폴더로 두세요 |
| `LICENSE` | 이용 조건 (비영리 학술·교육 목적 허용, 인용 필수) |
| `CITATION.cff` | 인용 정보 — 우측 **Cite this repository** 버튼에 연결됩니다 |

---

## 인용 · 권리

본 코드에 대한 권리 및 관련 특허는 **부경대학교 대기원격탐사연구실**
(Atmospheric Remote Sensing Laboratory, Pukyong National University) 에 있습니다.

본 분석 방법을 사용하시는 경우 아래 논문을 **반드시 인용**해 주시기 바랍니다.

> Shin, J., Kim, G., Kim, D., Tesche, M., Park, G., and Noh, Y.:
> Multi-section reference value for the analysis of horizontally scanning
> aerosol lidar observations,
> *Atmospheric Measurement Techniques*, **17**, 397–406, 2024.
> https://doi.org/10.5194/amt-17-397-2024

자세한 이용 조건은 [LICENSE](LICENSE) 를 참고하세요.

---

## 시작하기

### 1. 내려받기

녹색 **Code** 버튼 → **Download ZIP** → 압축 해제

### 2. 예제 데이터 준비

`testData.zip` 의 압축을 풀어, `.asc` 파일들이 **`testData`** 라는 이름의
폴더 안에 오도록 배치합니다. 코드가 그 폴더를 자동으로 찾습니다.

```
저장소 폴더/
├── Klett_Main_PKNU.m
├── Klett_Main_PKNU.py
└── testData/
    ├── A_2022-02-24-22-02-50_Lat-36.985476_..._Azimuth--42_Elevation-2.asc
    ├── B_2022-02-24-22-02-50_...
    └── ...
```

폴더 이름이 다르면 `testData` 로 바꾸거나, 코드 상단의 `DATA_DIR` 한 줄을
실제 폴더 이름으로 고치면 됩니다.

### 3-A. MATLAB 으로 실행

필요 환경

- MATLAB R2019b 이상 (`mean(...,'all')`, `sgtitle`, `writematrix` 사용)
- **Curve Fitting Toolbox** — 배경 신호의 `exp2` 피팅에 필요

```matlab
>> Klett_Main_PKNU
```

`testData` 폴더를 자동으로 찾아 분석 가능한 측정 시각을 목록으로 보여줍니다.

```
===== 분석 가능한 날짜 · 시각 =====
  1) 2022-02-24 22:02:50
  2) 2022-02-27 20:14:56
  3) 2022-03-09 20:14:01
```

코드 상단의 `SCAN_ID` 번호만 바꾸면 그 시각으로 전체 분석과 그림이 다시 돌아갑니다.
에디터에서 `Ctrl+Enter` 로 `%%` 구간을 하나씩 실행하면 단계별로 결과를
확인하며 진행할 수 있습니다.

### 3-B. Python 으로 실행

필요 환경 — `numpy`, `scipy`, `matplotlib`, `pandas`, `openpyxl`

```bash
pip install numpy scipy matplotlib pandas openpyxl
```

```bash
python Klett_Main_PKNU.py                  # 옆의 testData 폴더 사용
python Klett_Main_PKNU.py <폴더> [SCAN_ID]   # 폴더 · 스캔 번호 직접 지정
```

`# %%` 셀로 나뉘어 있어 VS Code · Spyder 에서 구간별 실행이 가능합니다.
MATLAB 판과 수치가 일치하도록 포팅했습니다
(열 우선 reshape, 음수 로그, NaN 전파 평균 등 동작을 맞춤).

### 4. 결과

`결과/` 폴더에 그림 7장(`01_raw.png` ~ `07_dp.png`)과
산출물(`products_*.mat`, `products_*.xlsx`)이 저장됩니다.

xlsx 는 시트마다 **1행 = 실제 방위각(deg), 1열 = 수평거리(km)** 형식입니다
(`EC_532P`, `EC_1064`, `AE`, `Dp`).

---

## 분석 흐름

```
0. 데이터 설정        경로 · 분석할 날짜/시각 · 채널
1. 데이터 불러오기     .asc -> [거리 x 방위각] 행렬      [그림 1] Raw 신호
2. 분석 고정값         해상도 · 참조거리 · k · 이득비
3. 거리 합산          6개씩 비중첩 합산 (4.8 -> 28.8 m, S/N √6배)
4. 배경 보정          두 추정값 중 작은 쪽 채택          [그림 2]
5. 이동평균 · RCS      S(r) = ln(P·r²)                  [그림 3]
6. 참조 소산계수       multi-section + 1σ 반복 제거      [그림 4]  ★ AMT 2024
7. Klett 소산계수      후방적분                          [그림 5]
8. Ångström exponent                                   [그림 6]
9. 편광비             Dp = S/(P+S)                      [그림 7]
10. 저장
```

### 이론 요약

**라이다 방정식** — 되돌아온 빛의 세기

$$P(r) = C \cdot O(r) \cdot \frac{\beta(r)}{r^2}\cdot\exp\left[-2\int_0^r \alpha(r')dr'\right]$$

미지수가 β(후방산란)와 α(소산) 둘인데 관측식은 하나뿐이라, 가정 없이는 풀 수 없습니다.

**거리보정신호** — r² 를 곱해 기하 감쇠를 지우고 ln 을 씌워 α 를 지수 밖으로 꺼냅니다.

$$S(r) \equiv \ln\left[P(r)r^2\right] = \ln\left[C\,O(r)\,\beta(r)\right] - 2\int_0^r\alpha(r')dr'$$

r 로 미분하면 시스템 상수 C 가 사라지고 적분이 풀립니다
(중첩이 완성된 구간, 즉 O = 1 에서):

$$\frac{dS}{dr} = \frac{1}{\beta}\frac{d\beta}{dr} - 2\alpha(r)$$

**두 가지 닫는 방법**

| | 가정 | 쓰이는 곳 |
|---|---|---|
| 기울기법 | β = const → α = −½·dS/dr | 참조 소산계수 α(rₘ) 산출 (6단계) |
| Klett | β = A·α^k → Bernoulli 형 | 전체 프로파일 (7단계) |

**Klett 해** — 참조점에서 라이다 쪽으로 거꾸로 적분하므로 경계조건 오차가 희석됩니다.

$$\alpha(r) = \frac{\exp\left[(S(r)-S(r_m))/k\right]}{\dfrac{1}{\alpha(r_m)} + \dfrac{2}{k}\displaystyle\int_{r}^{r_m}\exp\left[(S(r')-S(r_m))/k\right]dr'}$$

**Ångström exponent** — 파장의존성이 입자 크기를 말해줍니다.

$$\mathrm{AE} = -\frac{\ln(\alpha_{532}/\alpha_{1064})}{\ln(532/1064)}$$

**편광비** — 비구형 입자가 편광을 흐트러뜨리는 정도.

$$\mathrm{Dp} = \frac{S_{\mathrm{cal}}}{P + S_{\mathrm{cal}}},\qquad S_{\mathrm{cal}} = \frac{S}{G}$$

---

## 데이터 형식

`testData.zip` 에 들어있는 `.asc` 는 헤더 없이 한 줄에 한 bin 값(광자수)만
들어있는 텍스트입니다.

**파일명 규칙**

```
<채널>_<YYYY-MM-DD-HH-MM-SS>_Lat<위도>_Lon<경도>_Azimuth-<방위각>_Elevation-<고도각>.asc
```

구분자가 `-` 이므로 값이 음수면 하이픈이 두 번 나옵니다
(`Azimuth--15` → −15°, `Azimuth-5` → +5°).

| 채널 접두어 | 파장 |
|---|---|
| `A` | 1064 nm |
| `B` | 532 nm P (평행) |
| `C` | 532 nm S (수직) |

**예제 데이터** — 당진 (36.985476°N, 126.745613°E) 2022년 2–3월 야간 수평 스캔
3회분입니다. 1회 스캔은 방위각 22방향(−42 ~ 165°), 고도각 0~2°, 채널 3개 = 66개 파일.

---

## ★ 사용 전 반드시 확인할 세 가지

코드 상단에 모여 있습니다. **기본값은 가정이므로 실제 장비·교정값으로 바꿔야 합니다.**

| 변수 | 기본값 | 설명 |
|---|---|---|
| `res_raw` | `4.8` m | 원자료 1 bin 의 거리폭. `.asc` 에 들어있지 않습니다. 최대 샘플링 32 MHz 기준 4.8 m (Shin et al., 2024). **다른 설정으로 받은 자료라면 반드시 교체.** 틀리면 거리축 전체와 Klett 참조거리가 통째로 어긋납니다. |
| `az_north` | `61.5` deg | 실제 북쪽에 해당하는 파일명 azimuth 값. 당진 장비는 스캔 범위 가운데가 북쪽입니다. 설치 기록이 있으면 그 값으로. |
| `gain_ratio` | `1.0` | 편광 채널 이득비. **교정 전이므로 Dp 는 상대값입니다.** 절대 편광비로 해석하지 마세요. |

---

## 산출이 기대는 가정과 한계

논문에 쓰실 때 함께 보고해야 하는 항목입니다.

1. **k = 0.7 고정** — 라이다비를 사실상 고정한 것입니다. 다만 본 방법은 α(rₘ) 을
   기울기법으로 직접 측정하므로 k 민감도가 낮습니다(당진 2022-02-24 자료에서
   k 를 0.67–1.0 으로 바꿔도 EC532 변동 약 2.4 %). 라이다비 편차가 큰 조건에서는
   k 조정을 검토하시기 바랍니다.
2. **참조구간 균질 가정** — 후보 구간이 균질하다고 보고 기울기를 취합니다.
   창 안에서 β 가 10 % 만 기울어도 α 가 50 % 넘게 틀어질 수 있어, 여러 창·여러
   방향에서 뽑아 이상치를 거르는 것이 필수입니다. 본 방법의 불확도는 대부분
   이 단계에서 결정됩니다.
3. **α(rₘ) 하나를 전 방위각에 공통 적용** — 방향별 대기 상태가 크게 다르면
   편향이 생깁니다.
4. **AE 의 532 는 P 채널만** — 엄밀히는 P+S 가 전체 532 입니다. 이득비 교정 후
   합산하는 것이 맞습니다.
5. **중첩(overlap) 미보정** — 근거리를 배경 피팅에서 제외하는 것으로 대체합니다.
6. **분자(Rayleigh) 성분 미분리** — Klett 단일성분 해이므로 산출되는 α 는
   에어로졸 + 분자의 **총 소산계수**입니다. 표준대기 기준 532 nm 0.0117 /km,
   1064 nm 0.0007 /km 이며, 고도 변화가 있는 관측에서는 반드시 분리해야 합니다.

또 하나, 이산 적분 때문에 r = rₘ 에서 α 가 α(rₘ) 과 정확히 일치하지 않고
약 1 % 낮게 나옵니다(적분항이 0 이 아니라 한 bin 만큼 남기 때문). 거리가
가까워질수록 희석되어 1 km 에서는 0.4 % 수준입니다.

---

## 코드 구조

단일 파일입니다. 본문 아래에 지역 함수 13개가 모여 있습니다(두 버전 동일).

| 함수 | 역할 |
|---|---|
| `list_scans` | 폴더에서 측정 날짜·시각 목록 추출 |
| `read_asc` | `.asc` → [거리 × 방위각] 행렬 |
| `bin_sum` | 거리방향 비중첩 합산 |
| `subtract_bg` / `subtract_bg_raw` | 배경 추정·차감 (합산본 / 원해상도) |
| `fit_bg` | `exp2` 배경 피팅 |
| `movmean_kdh` | 양끝 비대칭 이동평균 |
| `ref_alpha` / `reference_alpha` | 기울기법 참조값 / multi-section 집계 |
| `reject_outlier` | 1σ 반복 이상치 제거 |
| `klett` | Klett 후방적분 |
| `show_product` | 산출물 3종 그림 (프로파일·단면·평면도) |
| `savepng` | 그림 저장 |

---

## 문의

부경대학교 대기원격탐사연구실 (Atmospheric Remote Sensing Laboratory,
Pukyong National University)
