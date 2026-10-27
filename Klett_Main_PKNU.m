%==========================================================================
%  Klett_Main_PKNU.m
%  스캐닝 라이다 수평 관측 — 소산계수 · AE · 편광비 산출 (단일 메인 코드)
%
%--------------------------------------------------------------------------
%  Copyright (C) 부경대학교 대기환경원격탐사연구실
%      Atmospheric Remote Sensing Laboratory,
%      Pukyong National University, Busan, Republic of Korea
%
%  본 코드에 대한 권리 및 관련 특허는 부경대학교 대기원격탐사연구실에 있습니다.
%
%  본 분석 방법(다중 구간 참조값 기반 수평 스캐닝 라이다 소산계수 산출)을
%  사용하시는 경우 아래 논문을 반드시 인용해 주시기 바랍니다.
%
%      Shin, J., Kim, G., Kim, D., Tesche, M., Park, G., and Noh, Y.:
%      Multi-section reference value for the analysis of horizontally
%      scanning aerosol lidar observations,
%      Atmospheric Measurement Techniques, 17, 397-406, 2024.
%      https://doi.org/10.5194/amt-17-397-2024
%--------------------------------------------------------------------------
%
%  흐름
%     0. 데이터 설정        경로 · 분석할 날짜/시각 · 채널
%     1. 데이터 불러오기     .asc -> [거리 x 방위각] 행렬
%        [그림 1] Raw 신호
%     2. 분석 고정값         여기만 고치면 전체가 따라감
%     3. 거리 합산          6개 비중첩 합산 (4.8 m -> 28.8 m)
%     4. 배경 보정          [그림 2]
%     5. 이동평균 · RCS      [그림 3]
%     6. 참조 소산계수       [그림 4]   <- AMT 2024 논문의 multi-section 방법
%     7. Klett 소산계수      [그림 5]
%     8. Angstrom exponent  [그림 6]
%     9. 편광비             [그림 7]
%    10. 저장
%

%
%  필요 툴박스 : Curve Fitting Toolbox (배경 피팅 exp2)
%==========================================================================
clear; clc; close all;

%% ===== 0. 데이터 설정 ====================================================
HERE     = fileparts(mfilename('fullpath'));
DATA_DIR = fullfile(HERE, 'testData');        % .asc 원자료 폴더
OUT_DIR  = fullfile(HERE, '결과');              % 그림 · mat 저장 폴더
if ~exist(OUT_DIR,'dir'); mkdir(OUT_DIR); end

% ---- 분석할 날짜 · 시각 고르기 ----------------------------------------
%  폴더 안의 스캔을 자동으로 찾아 목록으로 보여준다.
%  아래 SCAN_ID 번호만 바꾸면 그 날짜·시각으로 전체 분석과 그림이 다시 돌아간다.
SCAN_LIST = list_scans(DATA_DIR);

fprintf('===== 분석 가능한 날짜 · 시각 =====\n');
for i = 1:numel(SCAN_LIST)
    fprintf('  %d) %s\n', i, regexprep(SCAN_LIST{i}, ...
            '^(\d{4}-\d{2}-\d{2})-(\d{2})-(\d{2})-(\d{2})$', '$1 $2:$3:$4'));
end

SCAN_ID = 1;                       % <-- 여기 번호를 바꾸세요(현재 1~4)

if SCAN_ID < 1 || SCAN_ID > numel(SCAN_LIST)
    error('SCAN_ID 는 1 ~ %d 사이여야 합니다.', numel(SCAN_LIST));
end
SCAN = SCAN_LIST{SCAN_ID};
fprintf('  -> 선택: %s\n\n', SCAN);

% 채널 접두어 <-> 파장
CH.A = '1064';      % 1064 nm
CH.B = '532P';      % 532 nm 평행 (parallel)
CH.C = '532S';      % 532 nm 수직 (perpendicular)

% ★ 원자료 1 bin 의 거리폭 [m]. .asc 파일에는 없으므로 장비 설정값으로
%   확인해야 한다. 틀리면 거리축 전체와 Klett 참조거리가 통째로 어긋난다.
%   최대 샘플링 32 MHz -> 4.8 m (Shin et al., 2024, AMT 17, 397-406).
res_raw = 4.8;

SAVE_PNG = true;    % 그림을 png 로 저장할지

%% ===== 1. 데이터 불러오기 ================================================
fprintf('===== 1. 데이터 불러오기 =====\n');

[RAW.A, az, el, LAT, LON] = read_asc(DATA_DIR, SCAN, 'A');
RAW.B = read_asc(DATA_DIR, SCAN, 'B');
RAW.C = read_asc(DATA_DIR, SCAN, 'C');

nang = numel(az);                 % 방위각 개수
nbin = size(RAW.B,1);             % 원자료 bin 개수

fprintf('  스캔 %s\n', SCAN);
fprintf('  방위각 %d개 : %s deg\n', nang, strjoin(cellstr(compose('%g',az)), ', '));
fprintf('  고도각 %g ~ %g deg  (사실상 수평 관측)\n', min(el), max(el));
fprintf('  원자료 %d bin x %d 방위각,  라이다 위치 %.6fN %.6fE\n', ...
        nbin, nang, LAT, LON);

%% ----- [그림 1] Raw 신호 -------------------------------------------------
%  아직 아무 보정도 하지 않은 원자료. 분석 파라미터를 정하기 전에
%  신호가 어디까지 살아있는지 눈으로 먼저 확인한다.
r_raw = (1:nbin)' * res_raw;      % 거리축 [m]

figure('Name','[1] Raw 신호','Color','w');
subplot(1,2,1)
  plot(r_raw/1000, RAW.B); grid on
  xlim([0 6]);
  xlabel('수평거리 (km)'); ylabel('P(r)');
  title('원자료 P(r) — 532P, 전 방위각');
subplot(1,2,2)
  plot(r_raw/1000, real(log(RAW.B .* r_raw.^2))); grid on
  xlim([0 6]); ylim([15 25]);
  xlabel('수평거리 (km)'); ylabel('ln(P\cdotr^2)');
  title('거리보정신호 ln(P\cdotr^2)');
sgtitle(sprintf('[그림 1] Raw 신호 — %s', SCAN), 'Interpreter','none');
savepng(SAVE_PNG, OUT_DIR, '01_raw');

%% ===== 2. 분석 고정값 ====================================================
%  여기 값만 바꾸면 이후 전체가 따라간다. (res_raw 는 0번에서 설정)

avg     = 6;                 % 거리방향 비중첩 합산 개수
res     = res_raw * avg;     % [m] 합산 후 거리 해상도 (= 28.8 m)
ave_num = 3;                 % 이동평균 반폭 [points] (중앙 13점 평균)

n_tail  = 200;               % 배경 통계에 쓸 뒤쪽 bin 개수 (합산 후 기준)
n_tail_raw = 1200;           % 편광비용 배경 통계 bin 개수 (원해상도 기준)
fov_m   = 500;               % [m] 중첩(overlap) 미완성 구간 — 배경 피팅에서 제외

r_start = 750;               % [m] 참조거리 후보 시작점
r_width = 250;               % [m] 후보 간격 (= 기울기 적합 반폭)
ref_num = 20;                % 후보 개수 -> 1000 ~ 5750 m
far_ref = 5000;              % [m] Klett 후방적분 시작거리 rm
AE_K    = 0.7;               % Klett k 지수 (beta = A*alpha^k)
sigma   = 1;                 % 참조값 이상치 제거 기준 [sigma]

gain_ratio = 1.0;            % ★ 편광 채널 보정비율 — 교정 전이면 1 (상대값)

% 산출물 그림의 색 / 세로축 범위 ([] 로 두면 데이터에 맞춰 자동)
CLIM_EXT = [0 1];            % 소산계수 (1/km)   — 그림 5
CLIM_AE  = [0 3];            % Angstrom exponent — 그림 6
CLIM_DP  = [0 0.3];          % 편광비            — 그림 7
az_north   = 61.5;           % [deg] 실제 북쪽에 해당하는 파일명 azimuth 값

% 파생값
nsum      = floor(nbin/avg);            % 합산 후 bin 개수
dist      = (1:nsum)' * res;            % [m] 합산 후 거리축
Ref_Index = fix(far_ref / res);         % Klett 참조점 bin
az_true   = az - az_north;              % 실제 방위각 (북쪽 0도, 시계방향)
d_km      = dist(1:Ref_Index)/1000;     % [km] 산출물 거리축
sel       = d_km >= 0.5 & d_km <= 3;    % 통계용 대표 구간 0.5-3 km

fprintf('\n===== 2. 분석 고정값 =====\n');
fprintf('  거리 해상도 : %.2f m -> 합산 후 %.2f m (%d개씩)\n', res_raw, res, avg);
fprintf('  참조 후보   : %g ~ %g m, %d개\n', r_start+r_width, r_start+r_width*ref_num, ref_num);
fprintf('  Klett       : rm = %g m (bin %d),  k = %g\n', far_ref, Ref_Index, AE_K);
fprintf('  방위각 기준 : 장비 %g deg = 북쪽 -> 실제 %g ~ %g deg\n', ...
        az_north, min(az_true), max(az_true));

%% ===== 3. 거리 합산 ======================================================
%  광자 계수는 Poisson 통계를 따르므로 avg개를 합하면 상대오차가 sqrt(avg) 배 준다.
%  대가는 거리 해상도 (4.8 m -> 28.8 m).
fprintf('\n===== 3. 거리 합산 =====\n');

SUM.B = bin_sum(RAW.B, avg, nsum);
SUM.A = bin_sum(RAW.A, avg, nsum);
fprintf('  %d bin -> %d bin\n', nbin, nsum);

%% ===== 4. 배경 보정 ======================================================
%  하늘 배경 + 암전류를 뺀다. 두 가지로 추정해 "작은 쪽"을 택한다.
%    방법 1 : 맨 뒤 n_tail개 bin 의 평균 (2 sigma 이상치는 NaN 처리 후)
%    방법 2 : exp2 피팅 a*exp(b*r) + c 의 상수항 c
%  배경을 과대 추정하면 원거리 신호가 0 이하로 내려가 로그가 발산하므로
%  보수적으로 작은 값을 쓴다.
fprintf('\n===== 4. 배경 보정 =====\n');

j0 = ceil(nang/2);          % 그림에서 대표로 볼 방위각 (가운데)

[BG.B, bgv.B, bgfit.B, bgtail.B] = subtract_bg(SUM.B, dist, n_tail, fov_m);
[BG.A, bgv.A, bgfit.A, bgtail.A] = subtract_bg(SUM.A, dist, n_tail, fov_m);
fprintf('  532P 배경 %.1f ~ %.1f counts,  1064 배경 %.1f ~ %.1f counts\n', ...
        min(bgv.B), max(bgv.B), min(bgv.A), max(bgv.A));

%% ----- [그림 2] 배경 추정 -----------------------------------------------
figure('Name','[2] 배경 보정','Color','w');
subplot(1,2,1)
  plot(dist/1000, SUM.B(:,j0), 'Color',[.7 .7 .7]); hold on
  yline(bgtail.B(j0), 'r--','LineWidth',1.4);
  yline(bgfit.B(j0),  'b-.','LineWidth',1.4);
  yline(bgv.B(j0),    'k-', 'LineWidth',1.8);
  grid on; ylim([0 bgv.B(j0)*20]);      % 배경선이 보이도록 y축 확대
  xlabel('수평거리 (km)'); ylabel('P(r)');
  legend('합산 신호', ...
         sprintf('방법1 뒤 %d bin 평균 = %.1f', n_tail, bgtail.B(j0)), ...
         sprintf('방법2 exp2 상수항 = %.1f', bgfit.B(j0)), ...
         sprintf('채택 (작은 쪽) = %.1f', bgv.B(j0)), 'Location','southwest');
  title('배경 추정 — 두 방법 비교 (y축 확대)');
subplot(1,2,2)
  plot(dist/1000, BG.B(:,j0), 'b'); grid on
  xlim([0 8]); yline(0,'k:');
  xlabel('수평거리 (km)'); ylabel('P(r) - P_{bg}');
  title('배경 차감 후');
sgtitle(sprintf('[그림 2] 배경 보정 — 532P, 방위각 %g deg', az_true(j0)));
savepng(SAVE_PNG, OUT_DIR, '02_background');

%% ===== 5. 이동평균 · RCS =================================================
%  기울기를 읽을 수 있을 만큼만 매끄럽게 만든 뒤
%  S(r) = ln(P*r^2) 형태로 바꾼다. 여기서부터가 Klett 의 입력이다.
fprintf('\n===== 5. 이동평균 · RCS =====\n');

MM.B  = movmean_kdh(BG.B, ave_num);
MM.A  = movmean_kdh(BG.A, ave_num);
RCS.B = real(log(MM.B .* dist.^2));
RCS.A = real(log(MM.A .* dist.^2));
fprintf('  이동평균 반폭 %d (중앙 %d점 평균)\n', ave_num, 2*ave_num+1);

%% ----- [그림 3] RCS ------------------------------------------------------
figure('Name','[3] 이동평균 · RCS','Color','w');
subplot(1,2,1)
  plot(dist/1000, BG.B(:,j0), 'Color',[.7 .7 .7]); hold on
  plot(dist/1000, MM.B(:,j0), 'r','LineWidth',1.3); grid on
  xlim([0 6]); xlabel('수평거리 (km)'); ylabel('신호');
  legend('배경 차감', '이동평균 후','Location','northeast');
  title('이동평균');
subplot(1,2,2)
  plot(dist/1000, RCS.B); grid on
  xlim([0 6]); xlabel('수평거리 (km)'); ylabel('S(r) = ln(P\cdotr^2)');
  title('RCS — 전 방위각 (기울기가 곧 소산)');
sgtitle('[그림 3] 이동평균 · 거리보정신호');
savepng(SAVE_PNG, OUT_DIR, '03_rcs');

%% ===== 6. 참조 소산계수 ==================================================
%  수평 관측에서 대기가 균질하면  dS/dr = -2*alpha  이므로
%  기울기의 절반이 소산계수가 된다.
%  참조거리 후보 ref_num개 x 방위각 nang개를 모두 구한 뒤
%  1 sigma 이상치를 반복 제거하고 남은 값의 평균을 경계조건으로 쓴다.
%  ※ 채널마다 신호가 다르므로 참조값도 채널별로 따로 구한다.
fprintf('\n===== 6. 참조 소산계수 =====\n');

r_point = r_start + r_width*(1:ref_num);        % 후보 거리 [m]
ridx    = fix(r_point / res);                   % 후보 bin

[RefAlpha.B, ALPHA_REF, KEPT] = reference_alpha(MM.B, ridx, r_width, res, sigma);
RefAlpha.A                    = reference_alpha(MM.A, ridx, r_width, res, sigma);

fprintf('  후보 %d개 -> 이상치 제거 후 %d개 잔존 (532P 기준)\n', ...
        numel(ALPHA_REF), sum(isfinite(KEPT(:))));
fprintf('  RefAlpha  532P = %.4f /km,  1064 = %.4f /km\n', ...
        RefAlpha.B*1000, RefAlpha.A*1000);

%% ----- [그림 4] 참조값 결정 ---------------------------------------------
figure('Name','[4] 참조 소산계수','Color','w');
subplot(1,2,1)
  plot(r_point/1000, ALPHA_REF*1000, 'o-','MarkerSize',3,'LineWidth',.6); hold on
  yline(RefAlpha.B*1000,'k','LineWidth',2);
  grid on; xlabel('참조거리 후보 (km)'); ylabel('\alpha_{ref} (1/km)');
  title(sprintf('후보 %d개 x 방위각 %d개', ref_num, nang));
subplot(1,2,2)
  histogram(ALPHA_REF(:)*1000, 40, 'FaceColor',[.7 .7 .7]); hold on
  histogram(KEPT(isfinite(KEPT))*1000, 40, 'FaceColor',[.1 .2 .5]);
  xline(RefAlpha.B*1000,'r','LineWidth',2);
  grid on; xlabel('\alpha_{ref} (1/km)'); ylabel('개수');
  legend('전체', sprintf('%g\\sigma 제거 후',sigma), ...
         sprintf('평균 %.3f /km',RefAlpha.B*1000), 'Location','northeast');
  title('반복 이상치 제거');
sgtitle('[그림 4] 경계조건 \alpha(r_m) 결정 — 532P');
savepng(SAVE_PNG, OUT_DIR, '04_refalpha');

%% ===== 7. Klett 소산계수 =================================================
%  beta = A*alpha^k 를 가정해 미지수를 alpha 하나로 줄이고,
%  참조점 rm 에서 라이다 쪽으로 거꾸로 적분한다 (후방적분 -> 안정).
%
%            exp[(S(r) - S(rm))/k]
%  a(r) = ---------------------------------------------
%         1/a(rm) + (2/k) * integral_r^rm exp[(S-S(rm))/k] dr'
fprintf('\n===== 7. Klett 소산계수 =====\n');

EXT.B = klett(RCS.B, Ref_Index, res, AE_K, RefAlpha.B);
EXT.A = klett(RCS.A, Ref_Index, res, AE_K, RefAlpha.A);

fprintf('  532P 중앙값 %.4f /km,  1064 중앙값 %.4f /km\n', ...
        median(EXT.B(:),'omitnan')*1000, median(EXT.A(:),'omitnan')*1000);

%% ----- [그림 5] 소산계수 -------------------------------------------------
show_product(EXT.B*1000, d_km, az_true, '소산계수 (1/km)', ...
             sprintf('[그림 5] %s nm 소산계수 — %s', CH.B, SCAN), 'jet', CLIM_EXT);
savepng(SAVE_PNG, OUT_DIR, '05_ext532');

%% ===== 8. Angstrom exponent =============================================
%  alpha(lambda) ~ lambda^(-AE)  ->  AE = -ln(a532/a1064) / ln(532/1064)
%  주의: 여기서 532 는 P 채널만이다 (엄밀히는 P+S 가 전체 532).
fprintf('\n===== 8. Angstrom exponent =====\n');

aS = EXT.B;  aL = EXT.A;
bad = ~(aS>0) | ~(aL>0);                  % 로그를 취할 수 없는 점 제외
aS(bad) = NaN;  aL(bad) = NaN;
AE = -log(aS ./ aL) / log(532/1064);

fprintf('  제외(소산계수<=0) %d / %d 점\n', sum(bad(:)), numel(bad));
fprintf('  AE 중앙값 (0.5-3 km) = %.3f\n', median(AE(sel,:), 'all', 'omitnan'));

show_product(AE, d_km, az_true, 'Angstrom exponent', ...
             sprintf('[그림 6] AE (532/1064) — %s', SCAN), 'turbo', CLIM_AE);
savepng(SAVE_PNG, OUT_DIR, '06_ae');

%% ===== 9. 편광비 =========================================================
%  Dp = S_cal / (P + S_cal),   S_cal = S / gain_ratio
%  편광비는 원해상도로 배경보정한 뒤 구하는 것이 원래 방식이고,
%  소산계수·AE 와 같이 보려고 동일 격자(avg개 합산)로도 만든다.
fprintf('\n===== 9. 편광비 =====\n');

Praw = subtract_bg_raw(RAW.B, res_raw, n_tail_raw, n_tail, fov_m);
Sraw = subtract_bg_raw(RAW.C, res_raw, n_tail_raw, n_tail, fov_m);

Ps = bin_sum(Praw, avg, nsum);
Ss = bin_sum(Sraw, avg, nsum) / gain_ratio;
DP = Ss ./ (Ps + Ss);
DP(Ps<=0 | Ss<=0) = NaN;                  % 신호가 0 이하면 비가 무의미
DP(~isfinite(DP)) = NaN;
DP = DP(1:Ref_Index,:);

fprintf('  gain ratio = %g (1 이면 미교정 = 상대값)\n', gain_ratio);
fprintf('  Dp 중앙값 (0.5-3 km) = %.4f,  제외된 점 %d / %d\n', ...
        median(DP(sel,:),'all','omitnan'), sum(isnan(DP(:))), numel(DP));

show_product(DP, d_km, az_true, 'Dp = S/(P+S)', ...
             sprintf('[그림 7] 532 nm 편광비 (gain %g) — %s', gain_ratio, SCAN), 'turbo', CLIM_DP);
savepng(SAVE_PNG, OUT_DIR, '07_dp');

%% ===== 10. 저장 ===================================================
fprintf('\n===== 10. 저장 =====\n');

% mat 저장
PROD = struct('scan',SCAN, 'distance_km',d_km, 'azimuth',az_true, ...
              'lat',LAT, 'lon',LON, 'RefAlpha_532P',RefAlpha.B, 'RefAlpha_1064',RefAlpha.A, ...
              'EC_532P',EXT.B*1000, 'EC_1064',EXT.A*1000, 'AE',AE, 'Dp',DP);
save(fullfile(OUT_DIR, sprintf('products_%s.mat',SCAN)), 'PROD');

% xlsx 저장 — 시트마다 1행 = 방위각, 1열 = 수평거리
xls = fullfile(OUT_DIR, sprintf('products_%s.xlsx',SCAN));
if exist(xls,'file'); delete(xls); end
sheets = {'EC_532P',EXT.B*1000; 'EC_1064',EXT.A*1000; 'AE',AE; 'Dp',DP};
for s = 1:size(sheets,1)
    writematrix([NaN, az_true(:)'; d_km, sheets{s,2}], xls, 'Sheet', sheets{s,1});
end

fprintf('  저장 완료: %s\n', OUT_DIR);
fprintf('\n===== 전체 완료 =====\n');


%==========================================================================
%  로컬 함수
%==========================================================================

function tags = list_scans(data_dir)
%LIST_SCANS  폴더 안의 .asc 파일명에서 측정 날짜·시각을 뽑아 목록으로 돌려준다.
%  파일명 형식: <채널>_<YYYY-MM-DD-HH-MM-SS>_Lat..._Lon..._Azimuth..._Elevation....asc
    F = dir(fullfile(data_dir, '*.asc'));
    if isempty(F)
        error('list_scans:none', '.asc 파일이 없습니다: %s', data_dir);
    end
    t = regexp({F.name}, ...
        '^[A-Za-z]+_(\d{4}-\d{2}-\d{2}-\d{2}-\d{2}-\d{2})_Lat', 'tokens', 'once');
    t = t(~cellfun(@isempty, t));
    if isempty(t)
        error('list_scans:none', '파일명 규칙에 맞는 .asc 파일이 없습니다: %s', data_dir);
    end
    tags = unique(cellfun(@(c) c{1}, t, 'UniformOutput', false));
    tags = sort(tags);
end

% -------------------------------------------------------------------------
function [DATA, az, el, lat, lon] = read_asc(data_dir, scan, ch)
%READ_ASC  한 스캔 · 한 채널의 .asc 를 [거리 x 방위각] 행렬로 읽는다.
%  파일명: <채널>_<시각>_Lat<위도>_Lon<경도>_Azimuth-<방위각>_Elevation-<고도각>.asc
%  구분자가 '-' 라서 값이 음수면 하이픈이 두 번 나온다 (Azimuth--15 -> -15).
    pat = ['^' ch '_' scan '_Lat-(?<lat>-?[\d.]+)_Lon-(?<lon>-?[\d.]+)_' ...
           'Azimuth-(?<az>-?[\d.]+)_Elevation-(?<el>-?[\d.]+)\.asc$'];
    F = dir(fullfile(data_dir, [ch '_*.asc']));

    az = []; el = []; V = {}; lat = NaN; lon = NaN;
    for i = 1:numel(F)
        t = regexp(F(i).name, pat, 'names', 'once');
        if isempty(t); continue; end
        az(end+1) = str2double(t.az);   %#ok<AGROW>
        el(end+1) = str2double(t.el);   %#ok<AGROW>
        lat = str2double(t.lat);  lon = str2double(t.lon);
        V{end+1} = readmatrix(fullfile(data_dir, F(i).name), 'FileType','text'); %#ok<AGROW>
    end
    if isempty(az)
        error('read_asc:none', '파일을 찾지 못했습니다 (scan=%s, ch=%s)', scan, ch);
    end

    [az, o] = sort(az);  el = el(o);  V = V(o);        % 방위각 오름차순
    n = min(cellfun(@numel, V));                        % 공통 bin 개수
    DATA = zeros(n, numel(az));
    for j = 1:numel(az)
        DATA(:,j) = V{j}(1:n);
    end
end

% -------------------------------------------------------------------------
function S = bin_sum(DATA, avg, nsum)
%BIN_SUM  거리방향으로 avg개씩 비중첩 합산.
    S = reshape(DATA(1:avg*nsum, :), avg, nsum, []);
    S = reshape(sum(S, 1, 'omitnan'), nsum, []);
end

% -------------------------------------------------------------------------
function [OUT, bg, bg_fit, bg_tail] = subtract_bg(SIX, dist, n_tail, fov_m)
%SUBTRACT_BG  배경 추정 후 차감 (합산된 신호용).
%  방법1 뒤쪽 n_tail개 평균(2sigma 이상치 제거 후) 과
%  방법2 exp2 피팅 상수항 중 작은 값을 배경으로 쓴다.
    [n, nang] = size(SIX);
    i0  = n - n_tail;
    fov = find(dist == max(dist(dist < fov_m+1000)), 1, 'last');

    bg = zeros(1,nang); bg_fit = zeros(1,nang); bg_tail = zeros(1,nang);
    OUT = zeros(n, nang);
    for j = 1:nang
        % 2 sigma 이상치 -> NaN
        m = mean(SIX(i0:end,j),'omitnan');  s = std(SIX(i0:end,j),'omitnan');
        col = SIX(:,j);
        tail = i0:n;
        col(tail(col(tail) > m + 2*s)) = NaN;

        % 방법2 : exp2 피팅 상수항
        c = fit_bg(dist(fov:end), col(fov:end));
        bg_fit(j) = c;
        if ~(c > 0); c = 10000; end

        % 방법1 : 꼬리 평균
        nm = mean(col(i0:end),'omitnan');
        bg_tail(j) = nm;

        bg(j)    = min(c, nm);              % 작은 쪽 채택
        OUT(:,j) = col - bg(j);
    end
end

% -------------------------------------------------------------------------
function OUT = subtract_bg_raw(DATA, res_raw, n_stat, n_clip, fov_m)
%SUBTRACT_BG_RAW  편광비용 배경 차감 (원해상도, 거리합산 없음).
%  원해상도에서는 200개 bin 이 약 1 km 밖에 안 되므로 배경 통계는
%  뒤쪽 n_stat(=1200)개로 구하고, 이상치 제거는 뒤쪽 n_clip(=200)개에 적용한다.
    [n, nang] = size(DATA);
    dist  = (1:n)' * res_raw;
    istat = n - n_stat;
    iclip = n - n_clip;
    fov   = find(dist == max(dist(dist < fov_m)), 1, 'last');

    OUT = zeros(n, nang);
    for j = 1:nang
        m = mean(DATA(istat:end,j),'omitnan');  s = std(DATA(istat:end,j));
        col = DATA(:,j);
        clipr = iclip:n;
        col(clipr(col(clipr) > m + 2*s)) = NaN;

        c = fit_bg(dist(fov:end), col(fov:end));
        if ~(c > 0); c = 10000; end
        c = min(c, mean(col(iclip:end),'omitnan'));
        OUT(:,j) = col - c;
    end
end

% -------------------------------------------------------------------------
function c = fit_bg(t, y)
%FIT_BG  exp2 모델 a*exp(b*x) + c*exp(d*x) 로 피팅하고 상수항을 돌려준다.
%  d 의 상·하한을 모두 0 으로 묶어 실질적으로 a*exp(b*x) + c 가 되게 했으므로
%  세 번째 계수 c 가 배경값이다.  (Curve Fitting Toolbox 필요)

    SHOW_FIT = true;    % <-- 피팅 확인 그림. 이 줄만 지우면(또는 주석처리) 꺼진다.

    [xData, yData] = prepareCurveData(t, y);
    ft = fittype('exp2');
    opts = fitoptions('Method','NonlinearLeastSquares');
    opts.Display    = 'Off';
    opts.Lower      = [-Inf -Inf -Inf 0];
    opts.Upper      = [ Inf  Inf  Inf 0];
    opts.Normalize  = 'on';
    %opts.StartPoint = [1.24352873079792 -5.8973602901141 -64.5514470568424 0];
    opts.StartPoint = [5 -1 5000 0];
    f = fit(xData, yData, ft, opts);
    v = coeffvalues(f);
    c = v(3);

    % ---- 피팅 확인 그림 : 호출할 때마다 같은 창에 겹쳐 그린다 ------------
    %  c 가 0 이하로 나온 피팅은 호출부에서 채택하지 않고 고정값으로 대체되는데,
    %  그 경우에도 피팅 자체는 파란 점선으로 그려서 왜 버려졌는지 보이게 한다.
    %  ※ Normalize='on' 이라 계수 a·b 로 직접 식을 세우면 틀린다.
    %    feval(f,x) 가 정규화를 되돌려 평가해 주므로 반드시 이걸로 그린다.
    if exist('SHOW_FIT','var')
        h = findobj('Type','figure','Name','fit_bg');
        if isempty(h)
            h  = figure('Name','fit_bg','Color','w');
            ax = axes('Parent',h); hold(ax,'on'); grid(ax,'on');
            set(ax,'YScale','log');
            xlabel(ax,'수평거리 (km)'); ylabel(ax,'P(r)  (로그축)');
            title(ax,'fit\_bg — exp2 배경 피팅 (호출할 때마다 누적)');
            % 범례용 더미
            plot(ax,nan,nan,'.','Color',[.80 .80 .80],'MarkerSize',14, ...
                 'DisplayName','피팅에 쓴 신호');
            plot(ax,nan,nan,'-','Color',[.85 .20 .20],'LineWidth',1.2, ...
                 'DisplayName','exp2 피팅  (c > 0, 채택)');
            plot(ax,nan,nan,'--','Color',[.20 .40 .85],'LineWidth',1.2, ...
                 'DisplayName','exp2 피팅  (c \leq 0, 고정값으로 대체)');
            plot(ax,nan,nan,'o','MarkerFaceColor',[.85 .20 .20], ...
                 'MarkerEdgeColor','none','DisplayName','상수항 c');
            legend(ax,'Location','northeast');
            ylim([1 100000])
        else
            h  = h(1);
            ax = findobj(h,'Type','axes');  ax = ax(1);
        end

        xk = xData/1000;
        yp = yData;          yp(yp <= 0) = NaN;      % 로그축에 못 올리는 점은 생략
        yq = feval(f,xData); yq(yq <= 0) = NaN;

        plot(ax, xk, yp, '.', 'Color',[.80 .80 .80], 'MarkerSize',4, ...
             'HandleVisibility','off');
        if c > 0
            plot(ax, xk, yq, '-', 'Color',[.85 .20 .20], 'LineWidth',1.0, ...
                 'HandleVisibility','off');
            plot(ax, xk(end), c, 'o', 'MarkerSize',5, ...
                 'MarkerFaceColor',[.85 .20 .20], 'MarkerEdgeColor','none', ...
                 'HandleVisibility','off');
        else
            plot(ax, xk, yq, '--', 'Color',[.20 .40 .85], 'LineWidth',1.0, ...
                 'HandleVisibility','off');
        end
        drawnow limitrate
    end
end
% -------------------------------------------------------------------------
function OUT = movmean_kdh(X, n)
%MOVMEAN_KDH  이동평균. 양끝은 비대칭으로 처리한다 (원본 방식 유지).
%    앞쪽  k = 1..n           : mean(X(k : k+n))
%    중앙  k = n+1..end-n-1   : mean(X(k-n : k+n))
%    뒤쪽  k = end-n..end     : mean(X(k-n : k))
    [N, nang] = size(X);
    OUT = zeros(N, nang);
    for j = 1:nang
        v = X(:,j);
        for k = 1:n;            OUT(k,j) = mean(v(k : k+n));   end
        for k = n+1 : N-n-1;    OUT(k,j) = mean(v(k-n : k+n)); end
        for k = N-n : N;        OUT(k,j) = mean(v(k-n : k));   end
    end
end

% -------------------------------------------------------------------------
function a = ref_alpha(P, ridx, r_width, res)
%REF_ALPHA  참조거리 주변에서 ln(P*r^2) 의 기울기를 구해 |기울기|/2 를 돌려준다.
%  수평 + 균질 가정에서  dS/dr = -2*alpha  이므로 기울기의 절반이 소산계수다.
    iw = fix(r_width / res);
    k  = (ridx-iw) : (ridx+iw);
    X  = (k(:) * res);
    Y  = log(P(k(:)) .* X.^2);

    K    = [X, ones(numel(X),1)];
    coef = (K' * K) \ (K' * Y);        % 최소제곱 직선 적합
    a    = abs(coef(1)) / 2;
end

% -------------------------------------------------------------------------
function [RefAlpha, ALPHA_REF, KEPT] = reference_alpha(MM, ridx, r_width, res, sigma)
%REFERENCE_ALPHA  참조거리 후보 x 방위각 전체에서 기울기법으로 alpha 를 구하고,
%  이상치를 반복 제거한 뒤 남은 값의 평균을 경계조건 alpha(rm) 으로 돌려준다.
    nang = size(MM,2);  nref = numel(ridx);
    ALPHA_REF = zeros(nref, nang);
    for j = 1:nang
        for h = 1:nref
            ALPHA_REF(h,j) = ref_alpha(MM(:,j), ridx(h), r_width, res);
        end
    end
    KEPT     = reject_outlier(ALPHA_REF, sigma);
    RefAlpha = mean(KEPT(:), 'omitnan');
end

% -------------------------------------------------------------------------
function R = reject_outlier(R, sigma)
%REJECT_OUTLIER  mean +- sigma*std 를 벗어나는 값을 NaN 으로 반복 제거.
%  mean/std 를 열(방위각)별로 구하므로 limit 은 1 x nang 벡터이고,
%  '|' 로 묶은 if 문은 MATLAB 규칙상 "모든 방위각 기준에서 벗어날 때"만 참이 된다.
%  (원본 동작 유지)
    [nref, nang] = size(R);
    for it = 1:500
        m = mean(R,'omitnan');  s = std(R,'omitnan');
        hi = m + s*sigma;  lo = m - s*sigma;
        for j = 1:nang
            for p = 1:nref
                if (lo > R(p,j)) | (hi < R(p,j))   %#ok<OR2>
                    R(p,j) = NaN;
                end
            end
        end
    end
end

% -------------------------------------------------------------------------
function E = klett(RCS, Ref_Index, res, k, a_ref)
%KLETT  Klett(1981) 후방적분 해.
%              exp[(S(r) - S(rm))/k]
%  a(r) = ------------------------------------------------
%         1/a(rm) + (2/k) * sum_{r}^{rm} exp[(S-S(rm))/k] * dr
    nang = size(RCS,2);
    UP   = zeros(Ref_Index, nang);
    UPS  = zeros(Ref_Index, nang);
    for j = 1:nang
        UP(:,j) = exp( (RCS(1:Ref_Index,j) - RCS(Ref_Index,j)) / k );
        for i = 1:Ref_Index
            UPS(i,j) = sum(UP(i:Ref_Index,j), 'omitnan') * res;
        end
    end
    E = real( UP ./ (1/a_ref + (2/k)*UPS) );
end

% -------------------------------------------------------------------------
function show_product(V, d_km, az_true, lab, ttl, cmap, clim)
%SHOW_PRODUCT  산출물 하나를 프로파일 · 단면 · 평면도 세 장으로 보여준다.
%  단면은 측정한 방위각만 칸으로 배치하고(방위각 간격이 불균일하므로),
%  평면도는 방위각 사이를 보간해 부채꼴을 채운다.
    if nargin >= 7 && ~isempty(clim)
        cl = clim;                              % 지정된 범위 사용
    else
        v = sort(V(isfinite(V)));               % 비어 있으면 1~99 분위 자동
        if isempty(v)
            cl = [0 1];
        else
            cl = [v(max(1,round(0.01*numel(v)))), v(max(1,round(0.99*numel(v))))];
        end
        if ~(cl(2) > cl(1)); cl = [cl(1), cl(1)+eps]; end
    end
    nang = numel(az_true);
    if strcmp(cmap,'turbo') && isempty(which('turbo')); cmap = 'jet'; end  % 구버전 대비

    figure('Name', ttl, 'Color','w', 'Position',[80 80 1180 360]);

    subplot(1,3,1)
      plot(d_km, V, 'LineWidth',.8); grid on
      ylim(cl); xlabel('수평거리 (km)'); ylabel(lab);
      title('방위각별 수평 프로파일');

    subplot(1,3,2)
      imagesc(1:nang, d_km, V); set(gca,'YDir','normal');
      set(gca,'XTick',1:nang,'XTickLabel',cellstr(compose('%g',az_true)),'XTickLabelRotation',90);
      caxis(cl); c = colorbar; c.Label.String = lab; colormap(gca, cmap);
      xlabel('방위각 (deg)'); ylabel('수평거리 (km)');
      title('방위각 x 수평거리 단면');

    subplot(1,3,3)
      azf = linspace(min(az_true), max(az_true), 361);
      Z   = interp1(az_true(:), V', azf(:), 'linear')';
      [AZ, D] = meshgrid(azf, d_km);
      pcolor(D.*sind(AZ), D.*cosd(AZ), Z); shading interp
      caxis(cl); c = colorbar; c.Label.String = lab; colormap(gca, cmap);
      axis equal; grid on; hold on
      rmax = max(d_km);  th = linspace(0,2*pi,361);
      for rr = 1:floor(rmax); plot(rr*cos(th), rr*sin(th), 'k:','LineWidth',.4); end
      plot(0,0,'k^','MarkerFaceColor','k');
      text(0, rmax*1.12,'N','HorizontalAlignment','center','FontWeight','bold');
      text(rmax*1.12, 0,'E','HorizontalAlignment','left','FontWeight','bold');
      xlim([-1 1]*rmax*1.2); ylim([-1 1]*rmax*1.2);
      xlabel('동쪽 (km)'); ylabel('북쪽 (km)');
      title('평면도 (방위각 보간)');

    sgtitle(ttl, 'Interpreter','none');
end

% -------------------------------------------------------------------------
function savepng(flag, out_dir, name)
%SAVEPNG  현재 figure 를 png 로 저장.
    if flag
        saveas(gcf, fullfile(out_dir, [name '.png']));
    end
end
