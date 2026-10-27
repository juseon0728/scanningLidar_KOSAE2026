#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =========================================================================
#  Klett_Main_PKNU.py
#  스캐닝 라이다 수평 관측 — 소산계수 · AE · 편광비 산출 (단일 메인 코드)
#  Klett_Main_PKNU.m 의 Python 이식본. 계산 내용은 MATLAB 판과 동일하다.
#
# -------------------------------------------------------------------------
#  Copyright (C) 부경대학교 대기원격탐사연구실
#      Atmospheric Remote Sensing Laboratory,
#      Pukyong National University, Busan, Republic of Korea
#
#  본 코드에 대한 권리 및 관련 특허는 부경대학교 대기원격탐사연구실에 있습니다.
#
#  본 분석 방법(다중 구간 참조값 기반 수평 스캐닝 라이다 소산계수 산출)을
#  사용하시는 경우 아래 논문을 반드시 인용해 주시기 바랍니다.
#
#      Shin, J., Kim, G., Kim, D., Tesche, M., Park, G., and Noh, Y.:
#      Multi-section reference value for the analysis of horizontally
#      scanning aerosol lidar observations,
#      Atmospheric Measurement Techniques, 17, 397-406, 2024.
#      https://doi.org/10.5194/amt-17-397-2024
# -------------------------------------------------------------------------
#
#  흐름
#     0. 데이터 설정        경로 · 분석할 날짜/시각 · 채널
#     1. 데이터 불러오기     .asc -> [거리 x 방위각] 행렬
#        [그림 1] Raw 신호
#     2. 분석 고정값         여기만 고치면 전체가 따라감
#     3. 거리 합산          6개 비중첩 합산 (4.8 m -> 28.8 m)
#     4. 배경 보정          [그림 2]
#     5. 이동평균 · RCS      [그림 3]
#     6. 참조 소산계수       [그림 4]   <- AMT 2024 논문의 multi-section 방법
#     7. Klett 소산계수      [그림 5]
#     8. Angstrom exponent  [그림 6]
#     9. 편광비             [그림 7]
#    10. 저장
#
#  VS Code / Spyder 에서 '# %%' 구간을 Shift+Enter 로 하나씩 실행하면
#  MATLAB 의 Ctrl+Enter 처럼 단계별로 확인하며 진행할 수 있습니다.
#
#  필요 패키지 : numpy, scipy, matplotlib  (xlsx 저장에 pandas + openpyxl)
#      pip install numpy scipy matplotlib pandas openpyxl
#
#  실행
#      python Klett_Main_PKNU.py                 # 옆의 testData 폴더 사용
#      python Klett_Main_PKNU.py <폴더> [SCAN_ID]  # 폴더·스캔 번호 직접 지정
# =========================================================================
import os
import re
import sys
import glob
import warnings

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.io import savemat

warnings.filterwarnings('ignore')          # MATLAB 처럼 경고는 조용히


# -------------------------------------------------------------------------
def set_korean_font():
    """그림의 한글이 깨지지 않도록 시스템에 있는 한글 글꼴을 찾아 지정."""
    import matplotlib.font_manager as fm
    have = {f.name for f in fm.fontManager.ttflist}
    for name in ('Malgun Gothic', 'AppleGothic', 'NanumGothic', 'NanumBarunGothic',
                 'Noto Sans CJK KR', 'Noto Sans KR', 'Noto Sans CJK JP', 'IPAGothic'):
        if name in have:
            matplotlib.rcParams['font.family'] = name
            break
    matplotlib.rcParams['axes.unicode_minus'] = False   # 음수 기호 깨짐 방지


set_korean_font()


# =========================================================================
#  함수  (MATLAB 판의 로컬 함수와 1:1 대응)
# =========================================================================
def list_scans(data_dir):
    """폴더 안의 .asc 파일명에서 측정 날짜·시각을 뽑아 목록으로 돌려준다.

    파일명 형식: <채널>_<YYYY-MM-DD-HH-MM-SS>_Lat..._Lon..._Azimuth..._Elevation....asc
    """
    files = sorted(glob.glob(os.path.join(data_dir, '*.asc')))
    if not files:
        raise FileNotFoundError('.asc 파일이 없습니다: %s' % data_dir)
    pat = re.compile(r'^[A-Za-z]+_(\d{4}-\d{2}-\d{2}-\d{2}-\d{2}-\d{2})_Lat')
    tags = {m.group(1) for m in (pat.match(os.path.basename(f)) for f in files) if m}
    if not tags:
        raise FileNotFoundError('파일명 규칙에 맞는 .asc 파일이 없습니다: %s' % data_dir)
    return sorted(tags)


def read_asc(data_dir, scan, ch):
    """한 스캔 · 한 채널의 .asc 를 [거리 x 방위각] 행렬로 읽는다.

    파일명: <채널>_<시각>_Lat<위도>_Lon<경도>_Azimuth-<방위각>_Elevation-<고도각>.asc
    구분자가 '-' 라서 값이 음수면 하이픈이 두 번 나온다 (Azimuth--15 -> -15).
    """
    pat = re.compile(r'^' + ch + '_' + re.escape(scan) +
                     r'_Lat-(-?[\d.]+)_Lon-(-?[\d.]+)'
                     r'_Azimuth-(-?[\d.]+)_Elevation-(-?[\d.]+)\.asc$')
    az, el, V = [], [], []
    lat = lon = np.nan
    for f in sorted(glob.glob(os.path.join(data_dir, ch + '_*.asc'))):
        m = pat.match(os.path.basename(f))
        if m is None:
            continue
        lat, lon = float(m.group(1)), float(m.group(2))
        az.append(float(m.group(3)))
        el.append(float(m.group(4)))
        V.append(np.loadtxt(f).ravel())
    if not az:
        raise FileNotFoundError('파일을 찾지 못했습니다 (scan=%s, ch=%s)' % (scan, ch))

    o = np.argsort(az, kind='stable')                   # 방위각 오름차순
    az = np.array(az)[o]
    el = np.array(el)[o]
    V = [V[i] for i in o]
    n = min(len(v) for v in V)                          # 공통 bin 개수
    DATA = np.column_stack([v[:n] for v in V]).astype(float)
    return DATA, az, el, lat, lon


def bin_sum(DATA, avg, nsum):
    """거리방향으로 avg개씩 비중첩 합산."""
    X = DATA[:avg * nsum, :]
    # MATLAB 의 reshape 는 열 우선이므로 order='F'
    return np.nansum(X.reshape(avg, nsum, -1, order='F'), axis=0)


def subtract_bg(SIX, dist, n_tail, fov_m):
    """배경 추정 후 차감 (합산된 신호용).

    방법1 뒤쪽 n_tail개 평균(2sigma 이상치 제거 후) 과
    방법2 exp2 피팅 상수항 중 작은 값을 배경으로 쓴다.
    """
    n, nang = SIX.shape
    i0 = n - n_tail - 1                                  # MATLAB i0 = n-n_tail (1-based)
    fov = int(np.flatnonzero(dist < fov_m + 1000)[-1])   # 피팅 시작 bin (0-based)

    bg = np.zeros(nang)
    bg_fit = np.zeros(nang)
    bg_tail = np.zeros(nang)
    OUT = np.zeros((n, nang))
    for j in range(nang):
        # 2 sigma 이상치 -> NaN
        m = np.nanmean(SIX[i0:, j])
        s = np.nanstd(SIX[i0:, j], ddof=1)
        col = SIX[:, j].copy()
        tail = np.arange(i0, n)
        col[tail[col[tail] > m + 2 * s]] = np.nan

        # 방법2 : exp2 피팅 상수항
        c = fit_bg(dist[fov:], col[fov:])
        bg_fit[j] = c
        if not (c > 0):
            c = 10000.0

        # 방법1 : 꼬리 평균
        nm = np.nanmean(col[i0:])
        bg_tail[j] = nm

        bg[j] = min(c, nm)                               # 작은 쪽 채택
        OUT[:, j] = col - bg[j]
    return OUT, bg, bg_fit, bg_tail


def subtract_bg_raw(DATA, res_raw, n_stat, n_clip, fov_m):
    """편광비용 배경 차감 (원해상도, 거리합산 없음).

    원해상도에서는 200개 bin 이 약 1 km 밖에 안 되므로 배경 통계는
    뒤쪽 n_stat(=1200)개로 구하고, 이상치 제거는 뒤쪽 n_clip(=200)개에 적용한다.
    """
    n, nang = DATA.shape
    dist = np.arange(1, n + 1) * res_raw
    istat = n - n_stat - 1
    iclip = n - n_clip - 1
    fov = int(np.flatnonzero(dist < fov_m)[-1])

    OUT = np.zeros((n, nang))
    for j in range(nang):
        m = np.nanmean(DATA[istat:, j])
        s = np.std(DATA[istat:, j], ddof=1)
        col = DATA[:, j].copy()
        clipr = np.arange(iclip, n)
        col[clipr[col[clipr] > m + 2 * s]] = np.nan

        c = fit_bg(dist[fov:], col[fov:])
        if not (c > 0):
            c = 10000.0
        c = min(c, np.nanmean(col[iclip:]))
        OUT[:, j] = col - c
    return OUT


# ---- fit_bg 전용 전역 ------------------------------------------------------
SHOW_FIT = True       # <-- 피팅 확인 그림. 이 줄을 False 로 바꾸면 꺼진다.
_FIT_AX = [None]      # 누적해서 그릴 축을 담아두는 자리


def fit_bg(t, y):
    """exp2 모델 a*exp(b*x) + c*exp(d*x) 로 피팅하고 상수항을 돌려준다.

    d 를 0 으로 묶어 실질적으로 a*exp(b*x) + c 가 되므로 세 번째 계수 c 가 배경값이다.
    MATLAB 의 fitoptions Normalize='on' 과 같게 x 를 (x-mean)/std 로 정규화한 뒤
    피팅한다. 따라서 계수 a·b 로 원래 x 에 대한 식을 직접 세우면 틀린다.
    """
    t = np.asarray(t, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    ok = np.isfinite(t) & np.isfinite(y)                 # prepareCurveData 와 동일
    xData, yData = t[ok], y[ok]

    mu, sd = xData.mean(), xData.std(ddof=1)             # Normalize='on'
    xn = (xData - mu) / sd

    def model(x, a, b, c):
        return a * np.exp(b * x) + c

    try:
        p, _ = curve_fit(model, xn, yData,
                         p0=[5.0, -1.0, 5000.0], maxfev=40000)
    except Exception:
        return np.nan
    c = float(p[2])

    # ---- 피팅 확인 그림 : 호출할 때마다 같은 창에 겹쳐 그린다 --------------
    #  c 가 0 이하로 나온 피팅은 호출부에서 채택하지 않고 고정값으로 대체되는데,
    #  그 경우에도 피팅 자체는 파란 점선으로 그려서 왜 버려졌는지 보이게 한다.
    if SHOW_FIT:
        ax = _FIT_AX[0]
        if ax is None or not plt.fignum_exists(ax.figure.number):
            fig = plt.figure('fit_bg', facecolor='w')
            ax = fig.add_subplot(111)
            ax.set_yscale('log')
            ax.grid(True, alpha=.4)
            ax.set_xlabel('수평거리 (km)')
            ax.set_ylabel('P(r)  (로그축)')
            ax.set_title('fit_bg — exp2 배경 피팅 (호출할 때마다 누적)')
            ax.plot([], [], '.', color=(.80, .80, .80), ms=14, label='피팅에 쓴 신호')
            ax.plot([], [], '-', color=(.85, .20, .20), lw=1.2, label='exp2 피팅  (c > 0, 채택)')
            ax.plot([], [], '--', color=(.20, .40, .85), lw=1.2, label='exp2 피팅  (c ≤ 0, 고정값으로 대체)')
            ax.plot([], [], 'o', mfc=(.85, .20, .20), mec='none', label='상수항 c')
            ax.legend(loc='upper right', fontsize=8)
            ax.set_ylim(1, 100000)
            _FIT_AX[0] = ax

        xk = xData / 1000.0
        yp = np.where(yData > 0, yData, np.nan)          # 로그축에 못 올리는 점은 생략
        yq = model(xn, *p)
        yq = np.where(yq > 0, yq, np.nan)

        ax.plot(xk, yp, '.', color=(.80, .80, .80), ms=4)
        if c > 0:
            ax.plot(xk, yq, '-', color=(.85, .20, .20), lw=1.0)
            ax.plot(xk[-1], c, 'o', ms=5, mfc=(.85, .20, .20), mec='none')
        else:
            ax.plot(xk, yq, '--', color=(.20, .40, .85), lw=1.0)
    return c


def movmean_kdh(X, n):
    """이동평균. 양끝은 비대칭으로 처리한다 (원본 방식 유지).

    앞쪽  k = 1..n         : mean(X(k : k+n))
    중앙  k = n+1..N-n-1   : mean(X(k-n : k+n))
    뒤쪽  k = N-n..N       : mean(X(k-n : k))
    ※ MATLAB 의 mean 과 같게 NaN 을 전파시킨다 (omitnan 아님).
    """
    N, nang = X.shape
    OUT = np.zeros((N, nang))
    for j in range(nang):
        v = X[:, j]
        for k in range(0, n):                 # MATLAB k = 1..n
            OUT[k, j] = v[k:k + n + 1].mean()
        for k in range(n, N - n - 1):         # MATLAB k = n+1..N-n-1
            OUT[k, j] = v[k - n:k + n + 1].mean()
        for k in range(N - n - 1, N):         # MATLAB k = N-n..N
            OUT[k, j] = v[k - n:k + 1].mean()
    return OUT


def ref_alpha(P, ridx, r_width, res):
    """참조거리 주변에서 ln(P*r^2) 의 기울기를 구해 |기울기|/2 를 돌려준다.

    수평 + 균질 가정에서 dS/dr = -2*alpha 이므로 기울기의 절반이 소산계수다.
    ※ P 가 음수면 MATLAB 의 log 는 복소수를 돌려주고 마지막에 abs 로 크기를
      취하므로, 같은 결과가 나오도록 복소수로 계산한다.
    """
    iw = int(r_width // res)
    k = np.arange(ridx - iw, ridx + iw + 1)              # 0-based 창
    X = (k + 1) * res                                    # MATLAB k*res (1-based)
    Y = np.log(P[k].astype(complex) * X ** 2)

    K = np.column_stack([X, np.ones(len(X))])
    coef = np.linalg.lstsq(K.astype(complex), Y, rcond=None)[0]
    return abs(coef[0]) / 2


def reference_alpha(MM, ridx, r_width, res, sigma):
    """참조거리 후보 x 방위각 전체에서 기울기법으로 alpha 를 구하고,
    이상치를 반복 제거한 뒤 남은 값의 평균을 경계조건 alpha(rm) 으로 돌려준다."""
    nang = MM.shape[1]
    nref = len(ridx)
    ALPHA_REF = np.zeros((nref, nang))
    for j in range(nang):
        for h in range(nref):
            ALPHA_REF[h, j] = ref_alpha(MM[:, j], ridx[h], r_width, res)
    KEPT = reject_outlier(ALPHA_REF, sigma)
    return np.nanmean(KEPT), ALPHA_REF, KEPT


def reject_outlier(R, sigma):
    """mean +- sigma*std 를 벗어나는 값을 NaN 으로 반복 제거.

    mean/std 를 열(방위각)별로 구하므로 limit 은 길이 nang 벡터이고,
    MATLAB 의 '|' + if 조합은 "모든 방위각 기준에서 벗어날 때"만 참이 된다.
    (원본 동작 유지 — np.all 로 같은 뜻)
    """
    R = R.copy()
    nref, nang = R.shape
    for _ in range(500):
        m = np.nanmean(R, axis=0)
        s = np.nanstd(R, axis=0, ddof=1)
        hi = m + s * sigma
        lo = m - s * sigma
        for j in range(nang):
            for p in range(nref):
                if np.all((lo > R[p, j]) | (hi < R[p, j])):
                    R[p, j] = np.nan
    return R


def klett(RCS, Ref_Index, res, k, a_ref):
    """Klett(1981) 후방적분 해.

                  exp[(S(r) - S(rm))/k]
      a(r) = ------------------------------------------------
             1/a(rm) + (2/k) * sum_{r}^{rm} exp[(S-S(rm))/k] * dr
    """
    nang = RCS.shape[1]
    UP = np.zeros((Ref_Index, nang))
    UPS = np.zeros((Ref_Index, nang))
    for j in range(nang):
        UP[:, j] = np.exp((RCS[:Ref_Index, j] - RCS[Ref_Index - 1, j]) / k)
        for i in range(Ref_Index):
            UPS[i, j] = np.nansum(UP[i:Ref_Index, j]) * res
    return np.real(UP / (1.0 / a_ref + (2.0 / k) * UPS))


def show_product(V, d_km, az_true, lab, ttl, cmap, clim, out_dir=None, name=None):
    """산출물 하나를 프로파일 · 단면 · 평면도 세 장으로 보여준다.

    단면은 측정한 방위각만 칸으로 배치하고(방위각 간격이 불균일하므로),
    평면도는 방위각 사이를 보간해 부채꼴을 채운다.
    """
    if clim is not None and len(clim) == 2:
        cl = tuple(clim)
    else:
        v = np.sort(V[np.isfinite(V)])
        if v.size == 0:
            cl = (0.0, 1.0)
        else:
            lo = v[max(0, int(round(0.01 * v.size)) - 1)]
            hi = v[max(0, int(round(0.99 * v.size)) - 1)]
            cl = (lo, hi if hi > lo else lo + 1e-12)
    nang = len(az_true)
    if cmap == 'turbo' and 'turbo' not in plt.colormaps():
        cmap = 'jet'

    fig, axs = plt.subplots(1, 3, figsize=(13.1, 4.0), facecolor='w')
    fig.canvas.manager.set_window_title(ttl)

    ax = axs[0]
    ax.plot(d_km, V, lw=.8)
    ax.grid(True, alpha=.4); ax.set_ylim(cl)
    ax.set_xlabel('수평거리 (km)'); ax.set_ylabel(lab)
    ax.set_title('방위각별 수평 프로파일')

    ax = axs[1]
    im = ax.imshow(V, origin='lower', aspect='auto', cmap=cmap,
                   vmin=cl[0], vmax=cl[1],
                   extent=[0.5, nang + 0.5, d_km[0], d_km[-1]])
    ax.set_xticks(np.arange(1, nang + 1))
    ax.set_xticklabels(['%g' % a for a in az_true], rotation=90, fontsize=7)
    cb = fig.colorbar(im, ax=ax); cb.set_label(lab)
    ax.set_xlabel('방위각 (deg)'); ax.set_ylabel('수평거리 (km)')
    ax.set_title('방위각 x 수평거리 단면')

    ax = axs[2]
    azf = np.linspace(np.min(az_true), np.max(az_true), 361)
    Z = np.empty((len(d_km), azf.size))
    for i in range(len(d_km)):
        Z[i, :] = np.interp(azf, az_true, V[i, :])       # 방위각 선형 보간
    AZ, D = np.meshgrid(azf, d_km)
    pc = ax.pcolormesh(D * np.sin(np.deg2rad(AZ)), D * np.cos(np.deg2rad(AZ)), Z,
                       shading='gouraud', cmap=cmap, vmin=cl[0], vmax=cl[1])
    cb = fig.colorbar(pc, ax=ax); cb.set_label(lab)
    ax.set_aspect('equal'); ax.grid(True, alpha=.3)
    rmax = float(np.max(d_km))
    th = np.linspace(0, 2 * np.pi, 361)
    for rr in range(1, int(np.floor(rmax)) + 1):
        ax.plot(rr * np.cos(th), rr * np.sin(th), 'k:', lw=.4)
    ax.plot(0, 0, 'k^', mfc='k')
    ax.text(0, rmax * 1.12, 'N', ha='center', fontweight='bold')
    ax.text(rmax * 1.12, 0, 'E', ha='left', fontweight='bold')
    ax.set_xlim(-rmax * 1.2, rmax * 1.2); ax.set_ylim(-rmax * 1.2, rmax * 1.2)
    ax.set_xlabel('동쪽 (km)'); ax.set_ylabel('북쪽 (km)')
    ax.set_title('평면도 (방위각 보간)')

    fig.suptitle(ttl)
    fig.tight_layout()
    if out_dir and name:
        savepng(True, out_dir, name, fig)
    return fig


def savepng(flag, out_dir, name, fig=None):
    """현재 figure 를 png 로 저장."""
    if flag:
        (fig or plt.gcf()).savefig(os.path.join(out_dir, name + '.png'),
                                   dpi=150, bbox_inches='tight')


# =========================================================================
def main():
    # %% ===== 0. 데이터 설정 ==============================================
    HERE = os.path.dirname(os.path.abspath(__file__))
    DATA_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'testData')
    OUT_DIR = os.path.join(HERE, '결과')
    os.makedirs(OUT_DIR, exist_ok=True)

    # ---- 분석할 날짜 · 시각 고르기 --------------------------------------
    SCAN_LIST = list_scans(DATA_DIR)

    print('===== 분석 가능한 날짜 · 시각 =====')
    for i, tag in enumerate(SCAN_LIST, 1):
        print('  %d) %s' % (i, re.sub(r'^(\d{4}-\d{2}-\d{2})-(\d{2})-(\d{2})-(\d{2})$',
                                      r'\1 \2:\3:\4', tag)))

    SCAN_ID = 2                      # <-- 여기 번호를 바꾸세요
    if len(sys.argv) > 2:
        SCAN_ID = int(sys.argv[2])

    if SCAN_ID < 1 or SCAN_ID > len(SCAN_LIST):
        raise ValueError('SCAN_ID 는 1 ~ %d 사이여야 합니다.' % len(SCAN_LIST))
    SCAN = SCAN_LIST[SCAN_ID - 1]
    print('  -> 선택: %s\n' % SCAN)

    # 채널 접두어 <-> 파장
    CH = {'A': '1064',      # 1064 nm
          'B': '532P',      # 532 nm 평행 (parallel)
          'C': '532S'}      # 532 nm 수직 (perpendicular)

    # ★ 원자료 1 bin 의 거리폭 [m]. .asc 파일에는 없으므로 장비 설정값으로
    #   확인해야 한다. 틀리면 거리축 전체와 Klett 참조거리가 통째로 어긋난다.
    #   최대 샘플링 32 MHz -> 4.8 m (Shin et al., 2024, AMT 17, 397-406).
    res_raw = 4.8

    SAVE_PNG = True     # 그림을 png 로 저장할지

    # %% ===== 1. 데이터 불러오기 ==========================================
    print('===== 1. 데이터 불러오기 =====')

    RAW = {}
    RAW['A'], az, el, LAT, LON = read_asc(DATA_DIR, SCAN, 'A')
    RAW['B'], *_ = read_asc(DATA_DIR, SCAN, 'B')

    # 532S(C) 채널은 자료에 따라 없을 수 있다. 없으면 9단계 편광비만 건너뛴다.
    # (이 try 블록만 지우면 MATLAB 판과 완전히 같은 동작 — 없으면 오류로 멈춤)
    try:
        RAW['C'], *_ = read_asc(DATA_DIR, SCAN, 'C')
        HAS_C = True
    except FileNotFoundError:
        HAS_C = False
        print('  [알림] 532S(C) 채널 파일이 없어 편광비(9단계)는 건너뜁니다.')

    nang = len(az)                  # 방위각 개수
    nbin = RAW['B'].shape[0]        # 원자료 bin 개수

    print('  스캔 %s' % SCAN)
    print('  방위각 %d개 : %s deg' % (nang, ', '.join('%g' % a for a in az)))
    print('  고도각 %g ~ %g deg  (사실상 수평 관측)' % (el.min(), el.max()))
    print('  원자료 %d bin x %d 방위각,  라이다 위치 %.6fN %.6fE'
          % (nbin, nang, LAT, LON))

    # %% ----- [그림 1] Raw 신호 -------------------------------------------
    #  아직 아무 보정도 하지 않은 원자료. 분석 파라미터를 정하기 전에
    #  신호가 어디까지 살아있는지 눈으로 먼저 확인한다.
    r_raw = (np.arange(1, nbin + 1) * res_raw).reshape(-1, 1)   # 거리축 [m]

    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), facecolor='w')
    fig.canvas.manager.set_window_title('[1] Raw 신호')
    axs[0].plot(r_raw / 1000, RAW['B'], lw=.8)
    axs[0].grid(True, alpha=.4); axs[0].set_xlim(0, 6)
    axs[0].set_xlabel('수평거리 (km)'); axs[0].set_ylabel('P(r)')
    axs[0].set_title('원자료 P(r) — 532P, 전 방위각')
    axs[1].plot(r_raw / 1000, np.real(np.log((RAW['B'] * r_raw ** 2).astype(complex))), lw=.8)
    axs[1].grid(True, alpha=.4); axs[1].set_xlim(0, 6); axs[1].set_ylim(15, 25)
    axs[1].set_xlabel('수평거리 (km)'); axs[1].set_ylabel('ln(P·r²)')
    axs[1].set_title('거리보정신호 ln(P·r²)')
    fig.suptitle('[그림 1] Raw 신호 — %s' % SCAN)
    fig.tight_layout()
    savepng(SAVE_PNG, OUT_DIR, '01_raw', fig)

    # %% ===== 2. 분석 고정값 ==============================================
    #  여기 값만 바꾸면 이후 전체가 따라간다. (res_raw 는 0번에서 설정)

    avg = 6                      # 거리방향 비중첩 합산 개수
    res = res_raw * avg          # [m] 합산 후 거리 해상도 (= 28.8 m)
    ave_num = 3                  # 이동평균 반폭 [points] (중앙 7점 평균)

    n_tail = 200                 # 배경 통계에 쓸 뒤쪽 bin 개수 (합산 후 기준)
    n_tail_raw = 1200            # 편광비용 배경 통계 bin 개수 (원해상도 기준)
    fov_m = 500                  # [m] 중첩(overlap) 미완성 구간 — 배경 피팅에서 제외

    r_start = 750                # [m] 참조거리 후보 시작점
    r_width = 250                # [m] 후보 간격 (= 기울기 적합 반폭)
    ref_num = 20                 # 후보 개수 -> 1000 ~ 5750 m
    far_ref = 5000               # [m] Klett 후방적분 시작거리 rm
    AE_K = 0.7                   # Klett k 지수 (beta = A*alpha^k)
    sigma = 1                    # 참조값 이상치 제거 기준 [sigma]

    gain_ratio = 1.0             # ★ 편광 채널 보정비율 — 교정 전이면 1 (상대값)

    # 산출물 그림의 색 / 세로축 범위 (None 으로 두면 데이터에 맞춰 자동)
    CLIM_EXT = (0, 1)            # 소산계수 (1/km)   — 그림 5
    CLIM_AE = (0, 3)             # Angstrom exponent — 그림 6
    CLIM_DP = (0, 0.3)           # 편광비            — 그림 7
    az_north = 61.5              # [deg] 실제 북쪽에 해당하는 파일명 azimuth 값

    # 파생값
    nsum = nbin // avg                          # 합산 후 bin 개수
    dist = np.arange(1, nsum + 1) * res         # [m] 합산 후 거리축
    Ref_Index = int(far_ref // res)             # Klett 참조점 bin
    az_true = az - az_north                     # 실제 방위각 (북쪽 0도, 시계방향)
    d_km = dist[:Ref_Index] / 1000              # [km] 산출물 거리축
    sel = (d_km >= 0.5) & (d_km <= 3)           # 통계용 대표 구간 0.5-3 km

    print('\n===== 2. 분석 고정값 =====')
    print('  거리 해상도 : %.2f m -> 합산 후 %.2f m (%d개씩)' % (res_raw, res, avg))
    print('  참조 후보   : %g ~ %g m, %d개'
          % (r_start + r_width, r_start + r_width * ref_num, ref_num))
    print('  Klett       : rm = %g m (bin %d),  k = %g' % (far_ref, Ref_Index, AE_K))
    print('  방위각 기준 : 장비 %g deg = 북쪽 -> 실제 %g ~ %g deg'
          % (az_north, az_true.min(), az_true.max()))

    # %% ===== 3. 거리 합산 ================================================
    #  광자 계수는 Poisson 통계를 따르므로 avg개를 합하면 상대오차가 sqrt(avg) 배 준다.
    #  대가는 거리 해상도 (4.8 m -> 28.8 m).
    print('\n===== 3. 거리 합산 =====')
    SUM = {'B': bin_sum(RAW['B'], avg, nsum),
           'A': bin_sum(RAW['A'], avg, nsum)}
    print('  %d bin -> %d bin' % (nbin, nsum))

    # %% ===== 4. 배경 보정 ================================================
    #  하늘 배경 + 암전류를 뺀다. 두 가지로 추정해 "작은 쪽"을 택한다.
    #    방법 1 : 맨 뒤 n_tail개 bin 의 평균 (2 sigma 이상치는 NaN 처리 후)
    #    방법 2 : exp2 피팅 a*exp(b*r) + c 의 상수항 c
    #  배경을 과대 추정하면 원거리 신호가 0 이하로 내려가 로그가 발산하므로
    #  보수적으로 작은 값을 쓴다.
    print('\n===== 4. 배경 보정 =====')

    j0 = int(np.ceil(nang / 2)) - 1        # 그림에서 대표로 볼 방위각 (가운데)

    BG, bgv, bgfit, bgtail = {}, {}, {}, {}
    BG['B'], bgv['B'], bgfit['B'], bgtail['B'] = subtract_bg(SUM['B'], dist, n_tail, fov_m)
    BG['A'], bgv['A'], bgfit['A'], bgtail['A'] = subtract_bg(SUM['A'], dist, n_tail, fov_m)
    print('  532P 배경 %.1f ~ %.1f counts,  1064 배경 %.1f ~ %.1f counts'
          % (bgv['B'].min(), bgv['B'].max(), bgv['A'].min(), bgv['A'].max()))

    # %% ----- [그림 2] 배경 추정 ------------------------------------------
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), facecolor='w')
    fig.canvas.manager.set_window_title('[2] 배경 보정')
    ax = axs[0]
    ax.plot(dist / 1000, SUM['B'][:, j0], color=(.7, .7, .7))
    ax.axhline(bgtail['B'][j0], color='r', ls='--', lw=1.4,
               label='방법1 뒤 %d bin 평균 = %.1f' % (n_tail, bgtail['B'][j0]))
    ax.axhline(bgfit['B'][j0], color='b', ls='-.', lw=1.4,
               label='방법2 exp2 상수항 = %.1f' % bgfit['B'][j0])
    ax.axhline(bgv['B'][j0], color='k', ls='-', lw=1.8,
               label='채택 (작은 쪽) = %.1f' % bgv['B'][j0])
    ax.grid(True, alpha=.4); ax.set_ylim(0, bgv['B'][j0] * 20)   # 배경선이 보이도록 y축 확대
    ax.set_xlabel('수평거리 (km)'); ax.set_ylabel('P(r)')
    ax.legend(loc='lower left', fontsize=8)
    ax.set_title('배경 추정 — 두 방법 비교 (y축 확대)')
    ax = axs[1]
    ax.plot(dist / 1000, BG['B'][:, j0], 'b')
    ax.grid(True, alpha=.4); ax.set_xlim(0, 8); ax.axhline(0, color='k', ls=':')
    ax.set_xlabel('수평거리 (km)'); ax.set_ylabel('P(r) - P_bg')
    ax.set_title('배경 차감 후')
    fig.suptitle('[그림 2] 배경 보정 — 532P, 방위각 %g deg' % az_true[j0])
    fig.tight_layout()
    savepng(SAVE_PNG, OUT_DIR, '02_background', fig)

    # %% ===== 5. 이동평균 · RCS ===========================================
    #  기울기를 읽을 수 있을 만큼만 매끄럽게 만든 뒤
    #  S(r) = ln(P*r^2) 형태로 바꾼다. 여기서부터가 Klett 의 입력이다.
    print('\n===== 5. 이동평균 · RCS =====')
    MM = {'B': movmean_kdh(BG['B'], ave_num),
          'A': movmean_kdh(BG['A'], ave_num)}
    dcol = dist.reshape(-1, 1)
    RCS = {'B': np.real(np.log((MM['B'] * dcol ** 2).astype(complex))),
           'A': np.real(np.log((MM['A'] * dcol ** 2).astype(complex)))}
    print('  이동평균 반폭 %d (중앙 %d점 평균)' % (ave_num, 2 * ave_num + 1))

    # %% ----- [그림 3] RCS -------------------------------------------------
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), facecolor='w')
    fig.canvas.manager.set_window_title('[3] 이동평균 · RCS')
    axs[0].plot(dist / 1000, BG['B'][:, j0], color=(.7, .7, .7), label='배경 차감')
    axs[0].plot(dist / 1000, MM['B'][:, j0], 'r', lw=1.3, label='이동평균 후')
    axs[0].grid(True, alpha=.4); axs[0].set_xlim(0, 6)
    axs[0].set_xlabel('수평거리 (km)'); axs[0].set_ylabel('신호')
    axs[0].legend(loc='upper right', fontsize=8); axs[0].set_title('이동평균')
    axs[1].plot(dist / 1000, RCS['B'], lw=.8)
    axs[1].grid(True, alpha=.4); axs[1].set_xlim(0, 6)
    axs[1].set_xlabel('수평거리 (km)'); axs[1].set_ylabel('S(r) = ln(P·r²)')
    axs[1].set_title('RCS — 전 방위각 (기울기가 곧 소산)')
    fig.suptitle('[그림 3] 이동평균 · 거리보정신호')
    fig.tight_layout()
    savepng(SAVE_PNG, OUT_DIR, '03_rcs', fig)

    # %% ===== 6. 참조 소산계수 ============================================
    #  수평 관측에서 대기가 균질하면  dS/dr = -2*alpha  이므로
    #  기울기의 절반이 소산계수가 된다.
    #  참조거리 후보 ref_num개 x 방위각 nang개를 모두 구한 뒤
    #  1 sigma 이상치를 반복 제거하고 남은 값의 평균을 경계조건으로 쓴다.
    #  ※ 채널마다 신호가 다르므로 참조값도 채널별로 따로 구한다.
    print('\n===== 6. 참조 소산계수 =====')

    r_point = r_start + r_width * np.arange(1, ref_num + 1)     # 후보 거리 [m]
    ridx = (r_point // res).astype(int) - 1                     # 후보 bin (0-based)

    RefAlpha = {}
    RefAlpha['B'], ALPHA_REF, KEPT = reference_alpha(MM['B'], ridx, r_width, res, sigma)
    RefAlpha['A'], _, _ = reference_alpha(MM['A'], ridx, r_width, res, sigma)

    print('  후보 %d개 -> 이상치 제거 후 %d개 잔존 (532P 기준)'
          % (ALPHA_REF.size, int(np.isfinite(KEPT).sum())))
    print('  RefAlpha  532P = %.4f /km,  1064 = %.4f /km'
          % (RefAlpha['B'] * 1000, RefAlpha['A'] * 1000))

    # %% ----- [그림 4] 참조값 결정 -----------------------------------------
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), facecolor='w')
    fig.canvas.manager.set_window_title('[4] 참조 소산계수')
    axs[0].plot(r_point / 1000, ALPHA_REF * 1000, 'o-', ms=3, lw=.6)
    axs[0].axhline(RefAlpha['B'] * 1000, color='k', lw=2)
    axs[0].grid(True, alpha=.4)
    axs[0].set_xlabel('참조거리 후보 (km)'); axs[0].set_ylabel('α_ref (1/km)')
    axs[0].set_title('후보 %d개 x 방위각 %d개' % (ref_num, nang))
    axs[1].hist(ALPHA_REF.ravel() * 1000, bins=40, color=(.7, .7, .7), label='전체')
    axs[1].hist(KEPT[np.isfinite(KEPT)] * 1000, bins=40, color=(.1, .2, .5),
                label='%gσ 제거 후' % sigma)
    axs[1].axvline(RefAlpha['B'] * 1000, color='r', lw=2,
                   label='평균 %.3f /km' % (RefAlpha['B'] * 1000))
    axs[1].grid(True, alpha=.4)
    axs[1].set_xlabel('α_ref (1/km)'); axs[1].set_ylabel('개수')
    axs[1].legend(loc='upper right', fontsize=8); axs[1].set_title('반복 이상치 제거')
    fig.suptitle('[그림 4] 경계조건 α(r_m) 결정 — 532P')
    fig.tight_layout()
    savepng(SAVE_PNG, OUT_DIR, '04_refalpha', fig)

    # %% ===== 7. Klett 소산계수 ===========================================
    #  beta = A*alpha^k 를 가정해 미지수를 alpha 하나로 줄이고,
    #  참조점 rm 에서 라이다 쪽으로 거꾸로 적분한다 (후방적분 -> 안정).
    print('\n===== 7. Klett 소산계수 =====')
    EXT = {'B': klett(RCS['B'], Ref_Index, res, AE_K, RefAlpha['B']),
           'A': klett(RCS['A'], Ref_Index, res, AE_K, RefAlpha['A'])}
    print('  532P 중앙값 %.4f /km,  1064 중앙값 %.4f /km'
          % (np.nanmedian(EXT['B']) * 1000, np.nanmedian(EXT['A']) * 1000))

    # %% ----- [그림 5] 소산계수 --------------------------------------------
    show_product(EXT['B'] * 1000, d_km, az_true, '소산계수 (1/km)',
                 '[그림 5] %s nm 소산계수 — %s' % (CH['B'], SCAN), 'jet', CLIM_EXT,
                 OUT_DIR if SAVE_PNG else None, '05_ext532')

    # %% ===== 8. Angstrom exponent ========================================
    #  alpha(lambda) ~ lambda^(-AE)  ->  AE = -ln(a532/a1064) / ln(532/1064)
    #  주의: 여기서 532 는 P 채널만이다 (엄밀히는 P+S 가 전체 532).
    print('\n===== 8. Angstrom exponent =====')
    aS = EXT['B'].copy()
    aL = EXT['A'].copy()
    bad = ~(aS > 0) | ~(aL > 0)               # 로그를 취할 수 없는 점 제외
    aS[bad] = np.nan
    aL[bad] = np.nan
    AE = -np.log(aS / aL) / np.log(532 / 1064)

    print('  제외(소산계수<=0) %d / %d 점' % (int(bad.sum()), bad.size))
    print('  AE 중앙값 (0.5-3 km) = %.3f' % np.nanmedian(AE[sel, :]))

    show_product(AE, d_km, az_true, 'Angstrom exponent',
                 '[그림 6] AE (532/1064) — %s' % SCAN, 'turbo', CLIM_AE,
                 OUT_DIR if SAVE_PNG else None, '06_ae')

    # %% ===== 9. 편광비 ====================================================
    #  Dp = S_cal / (P + S_cal),   S_cal = S / gain_ratio
    #  편광비는 원해상도로 배경보정한 뒤 구하는 것이 원래 방식이고,
    #  소산계수·AE 와 같이 보려고 동일 격자(avg개 합산)로도 만든다.
    print('\n===== 9. 편광비 =====')
    if not HAS_C:
        DP = None
        print('  532S 채널이 없어 건너뜁니다.')
    else:
        Praw = subtract_bg_raw(RAW['B'], res_raw, n_tail_raw, n_tail, fov_m)
        Sraw = subtract_bg_raw(RAW['C'], res_raw, n_tail_raw, n_tail, fov_m)

        Ps = bin_sum(Praw, avg, nsum)
        Ss = bin_sum(Sraw, avg, nsum) / gain_ratio
        DP = Ss / (Ps + Ss)
        DP[(Ps <= 0) | (Ss <= 0)] = np.nan     # 신호가 0 이하면 비가 무의미
        DP[~np.isfinite(DP)] = np.nan
        DP = DP[:Ref_Index, :]

        print('  gain ratio = %g (1 이면 미교정 = 상대값)' % gain_ratio)
        print('  Dp 중앙값 (0.5-3 km) = %.4f,  제외된 점 %d / %d'
              % (np.nanmedian(DP[sel, :]), int(np.isnan(DP).sum()), DP.size))

        show_product(DP, d_km, az_true, 'Dp = S/(P+S)',
                     '[그림 7] 532 nm 편광비 (gain %g) — %s' % (gain_ratio, SCAN),
                     'turbo', CLIM_DP, OUT_DIR if SAVE_PNG else None, '07_dp')

    # %% ===== 10. 저장 =====================================================
    print('\n===== 10. 저장 =====')

    # mat 저장 (MATLAB 에서 load 해서 그대로 쓸 수 있음)
    PROD = {'scan': SCAN, 'distance_km': d_km.reshape(-1, 1),
            'azimuth': az_true.reshape(1, -1), 'lat': LAT, 'lon': LON,
            'RefAlpha_532P': RefAlpha['B'], 'RefAlpha_1064': RefAlpha['A'],
            'EC_532P': EXT['B'] * 1000, 'EC_1064': EXT['A'] * 1000, 'AE': AE}
    if HAS_C:
        PROD['Dp'] = DP
    savemat(os.path.join(OUT_DIR, 'products_%s.mat' % SCAN), {'PROD': PROD})

    # xlsx 저장 — 시트마다 1행 = 방위각, 1열 = 수평거리
    sheets = [('EC_532P', EXT['B'] * 1000), ('EC_1064', EXT['A'] * 1000), ('AE', AE)]
    if HAS_C:
        sheets.append(('Dp', DP))
    xls = os.path.join(OUT_DIR, 'products_%s.xlsx' % SCAN)
    try:
        import pandas as pd
        with pd.ExcelWriter(xls, engine='openpyxl') as xw:
            for nm, V in sheets:
                top = np.concatenate([[np.nan], az_true])
                body = np.column_stack([d_km, V])
                pd.DataFrame(np.vstack([top, body])).to_excel(
                    xw, sheet_name=nm, header=False, index=False)
    except ImportError:
        for nm, V in sheets:                      # pandas 가 없으면 csv 로
            top = np.concatenate([[np.nan], az_true])
            body = np.column_stack([d_km, V])
            np.savetxt(os.path.join(OUT_DIR, 'products_%s_%s.csv' % (SCAN, nm)),
                       np.vstack([top, body]), delimiter=',', fmt='%.6g')
        print('  pandas/openpyxl 이 없어 csv 로 저장했습니다.')

    print('  저장 완료: %s' % OUT_DIR)
    print('\n===== 전체 완료 =====')

    plt.show()


if __name__ == '__main__':
    main()
