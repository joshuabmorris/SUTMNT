


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% ╔═════════════════════════════════════════════════════════════════════════╗ %
% ║																			║ %
% ║    ┏━━━━━┳━━━┳━━━┓														║ %
% ║    ┗━┓ ┏━┫ ┃ ┃ ┃ ┃    Thermal-Mechanical Nonlinearity Testing (TMNT)    ║ %
% ║      ┃ ┃ ┃ ┃ ┃ ┃ ┃    Consolidated Experimental Post-processing Script  ║ %
% ║      ┣━┻━┫ ┣━━━┻━┫    (16/10/2026)										║ %
% ║      ┃ ┃ ┃ ┣━┓ ┏━┛														║ %
% ║      ┃ ┃ ┃ ┃ ┃ ┃	  Joshua Morris, 952911@swansea.ac.uk				║ %
% ║      ┗━┻━━━┛ ┗━┛	  Swansea University								║ %
% ║																			║ %
% ║    Supervisors:															║ %
% ║    Matthew Bonney,	m.s.bonney@swansea.ac.uk							║ %
% ║    Hamed Haddad Khodaparast, h.haddadkhodaparast@swansea.ac.uk			║ %
% ║																			║ %
% ║    MATLAB script with local functions: MATLAB R2016b or later.			║ %
% ║    The numerical analysis is intentionally kept in one script			║ %
% ║    with local helper functions so that the complete workflow			║ %
% ║    can be run from a single entry point.								║ %
% ║																			║ %
% ║    Required project layout:												║ %
% ║      data/                       Experimental CSV files					║ %
% ║																			║ %
% ║    Main outputs:														║ %
% ║      logs/                       Timestamped run log					║ %
% ║      exports/                    Numerical tables/intermediate exports	║ %
% ║      exports/validation/         Model-validation/robustness outputs	║ %
% ║      exports/journal/            Manuscript figures and tables			║ %
% ║																			║ %
% ╚═════════════════════════════════════════════════════════════════════════╝ %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Workflow overview:
%   1. Impact-hammer processing for monolithic and jointed specimens.
%   2. Stepped-sine processing over temperature and excitation level.
%   3. Forced-response resonance-peak extraction and displacement estimate.
%   4. Temperature-specific cubic backbone identification.
%   5. Thermal/model-comparison diagnostics and validation exports.
%   6. Standard-versus-extended empirical backbone comparison.
%   7. Journal figure/table generation from validated numerical exports.
%
% Interpretation notes used throughout the manuscript:
%   * Accelerometer temperature-sensitivity compensation was applied upstream
%     in LabVIEW. MATLAB receives already-corrected acceleration data.
%   * Stepped-sine resonance loci are forced-response peak loci, not direct
%     nonlinear-normal-mode backbones.
%   * Displacement amplitude is reconstructed from harmonic acceleration and
%     is not a directly measured displacement or modal coordinate.
%   * K* = m(2*pi*f)^2 uses physical specimen mass and is a stiffness-related
%     comparison only; it is not a formally identified modal stiffness.
%   * The A^4 backbone term is an empirical curvature term, not evidence of a
%     physical fifth-order restoring-force law.
%   * Residual-bootstrap intervals quantify conditional fitting sensitivity,
%     not full experimental uncertainty or repeatability.
%   * Energy-related quantities retained in the legacy processing are
%     uncalibrated proxies and are not headline dissipation measurements.


clc
clear variables
close all

if ~isfolder('data')
    error('research_pipeline:MissingData', ...
        'The data/ folder containing experimental CSV files is required.');
end

config.timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
if ~exist('logs','dir'),           mkdir('logs');           end
if ~exist('exports','dir'),        mkdir('exports');        end
if ~exist('exports/figures','dir'),mkdir('exports/figures');end

if ~exist(fullfile('exports','figures','paper'),'dir')
    mkdir(fullfile('exports','figures','paper'));
end

diary(fullfile('logs', ['log_' config.timestamp '.txt']));
fprintf('TMNT Post-processing  %s\n\n', config.timestamp);

config.steppedsine.nTemps = 16;
proc_all = struct();

% ── Step 1: Impact hammer ──────────────────────────────────────────────────
fprintf('Step 1: Impact hammer\n');

config.system.label    = 'Monolithic';
config.system.BeamType = 2;
config.system.Mass	   = 0.22821;
config.impact.data.fileName = 'H1xy_monolithic.csv';

%% ==== BEGIN ImpactHammer2.m (integrated) ====

config.impact.data.Folder = 'data';
config.impact.data.File   = fullfile(config.impact.data.Folder, config.impact.data.fileName);
config.impact.opts = detectImportOptions(config.impact.data.File);
config.impact.opts.DataLines = [14 Inf];

fprintf('Reading %s...\n', config.impact.data.fileName);
data = readmatrix(config.impact.data.File, config.impact.opts);
fprintf('...done.\n\n');

% Four-channel average (cols 2,4,6,8 are FRF magnitudes)
freq_Hz     = data(:,1);
chan_matrix = [data(:,2), data(:,4), data(:,6), data(:,8)];
valid_chans = ~all(chan_matrix == 0, 1) & ~all(isnan(chan_matrix), 1);
if sum(valid_chans) < 4
    warning('ImpactHammer2: %d of 4 channels valid.', sum(valid_chans));
end
resp_mean    = mean(chan_matrix(:, valid_chans), 2);
n_chans_used = sum(valid_chans);
fprintf('FRF averaged over %d of 4 channels.\n', n_chans_used);

% PCHIP sub-step resonance and half-power bandwidth
validF = ~isnan(resp_mean) & resp_mean > 0;
fD_imp = linspace(min(freq_Hz(validF)), max(freq_Hz(validF)), ...
                  max(5000, round((max(freq_Hz(validF))-min(freq_Hz(validF)))/0.001)));
rD_imp = ppval(pchip(freq_Hz(validF), resp_mean(validF)), fD_imp);

[peakAmp, iPk_imp] = max(rD_imp);
peakFreq  = fD_imp(iPk_imp);
halfPower = peakAmp / sqrt(2);

f1 = findHalfPowerImpact(fD_imp, rD_imp, iPk_imp, 'left',  halfPower);
f2 = findHalfPowerImpact(fD_imp, rD_imp, iPk_imp, 'right', halfPower);

if isnan(f1) || isnan(f2)
    warning('ImpactHammer2: could not find half-power crossings.');
    zeta_hpb = NaN;
else
    zeta_hpb = (f2 - f1) / (2*peakFreq);
end

fprintf('Peak (PCHIP): %.4f Hz\n', peakFreq);
% fprintf('HPB: f1=%.4f  f2=%.4f  zeta=%.5f\n', f1, f2, zeta_hpb);

% Nyquist circle fit (if phase data present in cols 3,5,7,9)
zeta_circle    = NaN;
circle_quality = NaN;

if size(data,2) >= 9
    try
        Hc_all = NaN(size(data,1), 4);
        for ch = 1:4
            if valid_chans(ch)
                Hc_all(:,ch) = data(:,2*ch) .* exp(1i * data(:,2*ch+1) * pi/180);
            end
        end
        Hc_mean  = mean(Hc_all(:,valid_chans), 2, 'omitnan');
        bw_est   = max(f2-f1, 4*mean(diff(freq_Hz)));
        win_mask = freq_Hz >= (peakFreq-3*bw_est) & freq_Hz <= (peakFreq+3*bw_est) & ...
                   ~isnan(Hc_mean) & abs(Hc_mean) > 0;
        Hc_win = Hc_mean(win_mask);
        f_win  = freq_Hz(win_mask);

        if numel(Hc_win) >= 5
            xp = real(Hc_win);  yp = imag(Hc_win);
            coef = [xp, yp, ones(size(xp))] \ -(xp.^2+yp.^2);
            xc = -coef(1)/2;  yc = -coef(2)/2;
            R2c = xc^2+yc^2-coef(3);
            if R2c > 0
                Rc = sqrt(R2c);
                circle_quality = mean(abs(sqrt((xp-xc).^2+(yp-yc).^2)-Rc)) / Rc;
                [~,i_near] = min(abs(f_win-peakFreq));
                i_lo = max(1,i_near-1);  i_hi = min(numel(f_win),i_near+1);
                if i_lo ~= i_hi && circle_quality < 0.3
                    dtheta = abs(angle(Hc_win(i_hi)-complex(xc,yc)) - angle(Hc_win(i_lo)-complex(xc,yc)));
                    if dtheta > pi; dtheta = 2*pi-dtheta; end
                    domega = 2*pi*(f_win(i_hi)-f_win(i_lo));
                    if domega > 0 && dtheta > 0
                        zc = 1 / (2*pi*peakFreq*(dtheta/domega));
                        if zc > 0.001 && zc < 0.15; zeta_circle = zc; end
                    end
                end
            end
        end
        % fprintf('Circle-fit: zeta=%.5f  quality=%.3f\n', zeta_circle, circle_quality);
    catch ME
        % fprintf('Circle fit skipped: %s\n', ME.message);
    end
else
    fprintf('Phase data not detected — circle fit skipped.\n');
end

% Best damping estimate
if ~isnan(zeta_circle)
    zeta = zeta_circle;  zeta_method = 'CircleFit';
else
    zeta = zeta_hpb;     zeta_method = 'HalfPower';
end
% fprintf('Damping: %.5f  (%s)\n\n', zeta, zeta_method);

% Modal mass — accelerance formula: m_eff = omega_n / (2*zeta*|H_peak|)
omega_n   = 2*pi*peakFreq;
m_eff_frf = 1 / (2*zeta*peakAmp);

if isfield(config,'steppedsine') && isfield(config.steppedsine,'system') && isfield(config.steppedsine.system,'Mass')
    m_phys = config.steppedsine.system.Mass;
elseif isfield(config,'system') && isfield(config.system,'Mass')
    m_phys = config.system.Mass;
else
    warning('ImpactHammer2: config.system.Mass not found — using FRF-derived m_eff.');
    m_phys = m_eff_frf;
end

% fprintf('m_eff FRF:  %.6f kg\n', m_eff_frf);
% fprintf('m_phys:     %.6f kg\n', m_phys);
if abs(m_eff_frf/m_phys - 1) > 0.25
    % warning('FRF m_eff differs from m_phys by %.1f%%.', 100*abs(m_eff_frf/m_phys-1));
end

% Physical-mass-based stiffness-related estimate:
%   m_eff_frf is retained as an FRF-derived diagnostic quantity. For the
%   cross-configuration comparison reported in the manuscript, the pipeline
%   deliberately uses the measured physical specimen mass in m*(2*pi*f_n)^2.
%   The resulting k_eff/K* quantity is therefore a stiffness-related scale,
%   not a formally identified modal stiffness. It should not be interpreted
%   as a rigorous modal parameter.
m_eff = m_phys;      % physical mass — consistent across systems for cross-comparison
k_eff = m_eff * omega_n^2;
Hpeak = peakAmp;
fprintf('k_eff (physical mass) = %.4f N/m  [m_phys=%.5f kg, m_eff_frf=%.5f kg]\n\n', k_eff, m_phys, m_eff_frf);

% Uncertainty
freq_step = mean(diff(freq_Hz));
sigma_f   = freq_step / 2;
sigma_k   = abs(2*m_eff*(2*pi)^2*peakFreq) * sigma_f;
zeta_uncert_frac = 0.10 + 0.10*strcmp(zeta_method,'HalfPower');
sigma_m_frf = m_eff_frf * zeta_uncert_frac;
% fprintf('Uncertainty:  sigma_f=+-%.4f Hz  sigma_k=+-%.1f N/m (%.2f%%)  sigma_m_frf=+-%.5f kg\n\n', ...
    % sigma_f, sigma_k, 100*sigma_k/k_eff, sigma_m_frf);

% Save to proc
result = struct('freq_Hz',freq_Hz,'response',resp_mean,'n_chans_used',n_chans_used, ...
    'peakFreq',peakFreq,'peakAmp',peakAmp,'Hpeak',Hpeak, ...
    'zeta',zeta,'zeta_hpb',zeta_hpb,'zeta_circle',zeta_circle, ...
    'zeta_method',zeta_method,'circle_quality',circle_quality, ...
    'm_eff',m_eff,'m_eff_frf',m_eff_frf,'m_phys',m_phys, ...
    'stiffness',k_eff, ...   % physical stiffness — comparable across systems
    'sigma_k',sigma_k,'sigma_k_pct',100*sigma_k/k_eff, ...
    'f1_hpb',f1,'f2_hpb',f2);

if strcmpi(config.impact.data.fileName, 'H1xy_monolithic.csv')
    proc.impact.monolithic = result;
    % fprintf('Monolithic:  f_n=%.4f Hz  zeta=%.5f (%s)  k=%.2f N/m\n\n', peakFreq,zeta,zeta_method,k_eff);
elseif strcmpi(config.impact.data.fileName, 'H1xy_jointed.csv')
    proc.impact.jointed = result;
    % fprintf('Jointed:     f_n=%.4f Hz  zeta=%.5f (%s)  k=%.2f N/m\n\n', peakFreq,zeta,zeta_method,k_eff);
else
    error('Unknown filename: expected H1xy_monolithic.csv or H1xy_jointed.csv');
end

% Verification plot (both systems overlaid on figure 1)
figure(1); hold on;
col = 'k'; if config.system.BeamType ~= 2; col = 'b'; end
plot(freq_Hz, resp_mean, [col '-'], 'LineWidth', 1.2);
plot(peakFreq, peakAmp, 'ro', 'MarkerFaceColor','r', 'MarkerSize', 8);
% if ~isnan(f1) && ~isnan(f2)
%     yline(halfPower,'m--','LabelHorizontalAlignment','left');
%     xline(f1,'m:','LineWidth',0.8);  xline(f2,'m:','LineWidth',0.8);
% end
xlabel('Frequency (Hz)'); ylabel('|H|  (m/s^2)/N');
grid on; set(gca,'YScale','log');

fprintf('Impact hammer processing complete.\n\n');

% =========================================================================
% LOCAL FUNCTION
% =========================================================================
%% ==== END ImpactHammer2.m ====

proc_all.Monolithic.impact = proc.impact;

config.system.label    = 'Jointed';
config.system.BeamType = 3;
config.system.Mass	   = 0.23304;
config.impact.data.fileName = 'H1xy_jointed.csv';

%% ==== BEGIN ImpactHammer2.m (integrated) ====

config.impact.data.Folder = 'data';
config.impact.data.File   = fullfile(config.impact.data.Folder, config.impact.data.fileName);
config.impact.opts = detectImportOptions(config.impact.data.File);
config.impact.opts.DataLines = [14 Inf];

fprintf('Reading %s...\n', config.impact.data.fileName);
data = readmatrix(config.impact.data.File, config.impact.opts);
fprintf('...done.\n\n');

% Four-channel average (cols 2,4,6,8 are FRF magnitudes)
freq_Hz     = data(:,1);
chan_matrix = [data(:,2), data(:,4), data(:,6), data(:,8)];
valid_chans = ~all(chan_matrix == 0, 1) & ~all(isnan(chan_matrix), 1);
if sum(valid_chans) < 4
    warning('ImpactHammer2: %d of 4 channels valid.', sum(valid_chans));
end
resp_mean    = mean(chan_matrix(:, valid_chans), 2);
n_chans_used = sum(valid_chans);
fprintf('FRF averaged over %d of 4 channels.\n', n_chans_used);

% PCHIP sub-step resonance and half-power bandwidth
validF = ~isnan(resp_mean) & resp_mean > 0;
fD_imp = linspace(min(freq_Hz(validF)), max(freq_Hz(validF)), ...
                  max(5000, round((max(freq_Hz(validF))-min(freq_Hz(validF)))/0.001)));
rD_imp = ppval(pchip(freq_Hz(validF), resp_mean(validF)), fD_imp);

[peakAmp, iPk_imp] = max(rD_imp);
peakFreq  = fD_imp(iPk_imp);
halfPower = peakAmp / sqrt(2);

f1 = findHalfPowerImpact(fD_imp, rD_imp, iPk_imp, 'left',  halfPower);
f2 = findHalfPowerImpact(fD_imp, rD_imp, iPk_imp, 'right', halfPower);

if isnan(f1) || isnan(f2)
    warning('ImpactHammer2: could not find half-power crossings.');
    zeta_hpb = NaN;
else
    zeta_hpb = (f2 - f1) / (2*peakFreq);
end

fprintf('Peak (PCHIP): %.4f Hz\n', peakFreq);
% fprintf('HPB: f1=%.4f  f2=%.4f  zeta=%.5f\n', f1, f2, zeta_hpb);

% Nyquist circle fit (if phase data present in cols 3,5,7,9)
zeta_circle    = NaN;
circle_quality = NaN;

if size(data,2) >= 9
    try
        Hc_all = NaN(size(data,1), 4);
        for ch = 1:4
            if valid_chans(ch)
                Hc_all(:,ch) = data(:,2*ch) .* exp(1i * data(:,2*ch+1) * pi/180);
            end
        end
        Hc_mean  = mean(Hc_all(:,valid_chans), 2, 'omitnan');
        bw_est   = max(f2-f1, 4*mean(diff(freq_Hz)));
        win_mask = freq_Hz >= (peakFreq-3*bw_est) & freq_Hz <= (peakFreq+3*bw_est) & ...
                   ~isnan(Hc_mean) & abs(Hc_mean) > 0;
        Hc_win = Hc_mean(win_mask);
        f_win  = freq_Hz(win_mask);

        if numel(Hc_win) >= 5
            xp = real(Hc_win);  yp = imag(Hc_win);
            coef = [xp, yp, ones(size(xp))] \ -(xp.^2+yp.^2);
            xc = -coef(1)/2;  yc = -coef(2)/2;
            R2c = xc^2+yc^2-coef(3);
            if R2c > 0
                Rc = sqrt(R2c);
                circle_quality = mean(abs(sqrt((xp-xc).^2+(yp-yc).^2)-Rc)) / Rc;
                [~,i_near] = min(abs(f_win-peakFreq));
                i_lo = max(1,i_near-1);  i_hi = min(numel(f_win),i_near+1);
                if i_lo ~= i_hi && circle_quality < 0.3
                    dtheta = abs(angle(Hc_win(i_hi)-complex(xc,yc)) - angle(Hc_win(i_lo)-complex(xc,yc)));
                    if dtheta > pi; dtheta = 2*pi-dtheta; end
                    domega = 2*pi*(f_win(i_hi)-f_win(i_lo));
                    if domega > 0 && dtheta > 0
                        zc = 1 / (2*pi*peakFreq*(dtheta/domega));
                        if zc > 0.001 && zc < 0.15; zeta_circle = zc; end
                    end
                end
            end
        end
        % fprintf('Circle-fit: zeta=%.5f  quality=%.3f\n', zeta_circle, circle_quality);
    catch ME
        % fprintf('Circle fit skipped: %s\n', ME.message);
    end
else
    fprintf('Phase data not detected — circle fit skipped.\n');
end

% Best damping estimate
if ~isnan(zeta_circle)
    zeta = zeta_circle;  zeta_method = 'CircleFit';
else
    zeta = zeta_hpb;     zeta_method = 'HalfPower';
end
% fprintf('Damping: %.5f  (%s)\n\n', zeta, zeta_method);

% Modal mass — accelerance formula: m_eff = omega_n / (2*zeta*|H_peak|)
omega_n   = 2*pi*peakFreq;
m_eff_frf = 1 / (2*zeta*peakAmp);

if isfield(config,'steppedsine') && isfield(config.steppedsine,'system') && isfield(config.steppedsine.system,'Mass')
    m_phys = config.steppedsine.system.Mass;
elseif isfield(config,'system') && isfield(config.system,'Mass')
    m_phys = config.system.Mass;
else
    warning('ImpactHammer2: config.system.Mass not found — using FRF-derived m_eff.');
    m_phys = m_eff_frf;
end

% fprintf('m_eff FRF:  %.6f kg\n', m_eff_frf);
% fprintf('m_phys:     %.6f kg\n', m_phys);
if abs(m_eff_frf/m_phys - 1) > 0.25
    % warning('FRF m_eff differs from m_phys by %.1f%%.', 100*abs(m_eff_frf/m_phys-1));
end

% Physical-mass-based stiffness-related estimate:
%   m_eff_frf is retained as an FRF-derived diagnostic quantity. For the
%   cross-configuration comparison reported in the manuscript, the pipeline
%   deliberately uses the measured physical specimen mass in m*(2*pi*f_n)^2.
%   The resulting k_eff/K* quantity is therefore a stiffness-related scale,
%   not a formally identified modal stiffness. It should not be interpreted
%   as a rigorous modal parameter.
m_eff = m_phys;      % physical mass — consistent across systems for cross-comparison
k_eff = m_eff * omega_n^2;
Hpeak = peakAmp;
fprintf('k_eff (physical mass) = %.4f N/m  [m_phys=%.5f kg, m_eff_frf=%.5f kg]\n\n', k_eff, m_phys, m_eff_frf);

% Uncertainty
freq_step = mean(diff(freq_Hz));
sigma_f   = freq_step / 2;
sigma_k   = abs(2*m_eff*(2*pi)^2*peakFreq) * sigma_f;
zeta_uncert_frac = 0.10 + 0.10*strcmp(zeta_method,'HalfPower');
sigma_m_frf = m_eff_frf * zeta_uncert_frac;
% fprintf('Uncertainty:  sigma_f=+-%.4f Hz  sigma_k=+-%.1f N/m (%.2f%%)  sigma_m_frf=+-%.5f kg\n\n', ...
    % sigma_f, sigma_k, 100*sigma_k/k_eff, sigma_m_frf);

% Save to proc
result = struct('freq_Hz',freq_Hz,'response',resp_mean,'n_chans_used',n_chans_used, ...
    'peakFreq',peakFreq,'peakAmp',peakAmp,'Hpeak',Hpeak, ...
    'zeta',zeta,'zeta_hpb',zeta_hpb,'zeta_circle',zeta_circle, ...
    'zeta_method',zeta_method,'circle_quality',circle_quality, ...
    'm_eff',m_eff,'m_eff_frf',m_eff_frf,'m_phys',m_phys, ...
    'stiffness',k_eff, ...   % physical stiffness — comparable across systems
    'sigma_k',sigma_k,'sigma_k_pct',100*sigma_k/k_eff, ...
    'f1_hpb',f1,'f2_hpb',f2);

if strcmpi(config.impact.data.fileName, 'H1xy_monolithic.csv')
    proc.impact.monolithic = result;
    % fprintf('Monolithic:  f_n=%.4f Hz  zeta=%.5f (%s)  k=%.2f N/m\n\n', peakFreq,zeta,zeta_method,k_eff);
elseif strcmpi(config.impact.data.fileName, 'H1xy_jointed.csv')
    proc.impact.jointed = result;
    % fprintf('Jointed:     f_n=%.4f Hz  zeta=%.5f (%s)  k=%.2f N/m\n\n', peakFreq,zeta,zeta_method,k_eff);
else
    error('Unknown filename: expected H1xy_monolithic.csv or H1xy_jointed.csv');
end

% Verification plot (both systems overlaid on figure 1)
figure(1); hold on;
col = 'k'; if config.system.BeamType ~= 2; col = 'b'; end
plot(freq_Hz, resp_mean, [col '-'], 'LineWidth', 1.2);
plot(peakFreq, peakAmp, 'ro', 'MarkerFaceColor','r', 'MarkerSize', 8);
% if ~isnan(f1) && ~isnan(f2)
%     yline(halfPower,'m--','LabelHorizontalAlignment','left');
%     xline(f1,'m:','LineWidth',0.8);  xline(f2,'m:','LineWidth',0.8);
% end
xlabel('Frequency (Hz)'); ylabel('|H|  (m/s^2)/N');
grid on; set(gca,'YScale','log');

fprintf('Impact hammer processing complete.\n\n');

% =========================================================================
% LOCAL FUNCTION
% =========================================================================
%% ==== END ImpactHammer2.m ====

proc_all.Jointed.impact = proc.impact;

% ── Step 2: Stepped-sine ──────────────────────────────────────────────────
fprintf('\nStep 2: Stepped-sine\n');

config.system.label    = 'Monolithic';
config.system.BeamType = 2;
config.steppedsine.data.fileName = arrayfun(@(k) ...
    sprintf("0000-00%02d_147-152-0.1_1000_2000_2.csv", k), 1:config.steppedsine.nTemps);

%% ==== BEGIN SteppedSine2.m (integrated) ====

clear c1_fit  % Preserve exportCache across systems for combined results.
config.steppedsine.temperature.Vector = -30:10:120;
proc.steppedsine.backbone = cell(config.steppedsine.nTemps, 1);
proc.steppedsine.amp      = cell(config.steppedsine.nTemps, 1);
proc.steppedsine.FRF      = cell(config.steppedsine.nTemps, 1);

fJ0 = proc.impact.jointed.peakFreq;
fM0 = proc.impact.monolithic.peakFreq;

% Pre-loop preallocation
nTemps_pre = numel(-30:10:120);
Na_pre     = numel(0.5:0.5:8);

AMatrix        = NaN(Na_pre, nTemps_pre);
kEffMatrix     = NaN(Na_pre, nTemps_pre);
cEqMatrix      = NaN(Na_pre, nTemps_pre);
zetaMatrix     = NaN(Na_pre, nTemps_pre);
peakFreqMatrix = NaN(Na_pre, nTemps_pre);
peakAmpMatrix  = NaN(Na_pre, nTemps_pre);

summaryNonlinearityType = cell(nTemps_pre, 1);
summaryOnset            = NaN(1, nTemps_pre);
summaryOnsetUm          = NaN(1, nTemps_pre);

% Legacy descriptive onset threshold
isMonolithic_early = strcmp(config.system.label, 'Monolithic');
if isMonolithic_early
    onset_thresh_pct = 0.20;   % legacy descriptive threshold for monolithic case
else
    onset_thresh_pct = 0.50;   % legacy descriptive threshold for jointed case
end

% =========================================================================
% TEMPERATURE LOOP
% =========================================================================

for Z = 1:config.steppedsine.nTemps

    config.steppedsine.data.Name = config.steppedsine.data.fileName(Z);
    config.steppedsine.data.File = fullfile('data', config.steppedsine.data.Name);
    fprintf('Reading %s...\n', config.steppedsine.data.Name);
    config.steppedsine.data.opts = detectImportOptions(config.steppedsine.data.File);
    config.steppedsine.data.opts.DataLines = [24, Inf];
    config.steppedsine.data.Table = readtable(config.steppedsine.data.File, config.steppedsine.data.opts);

    % Parse metadata from filename
    p = sscanf(config.steppedsine.data.Name, '%d-%d_%d-%d-%f_%d_%d_%d');
    if numel(p) < 8; error('Filename parsing failed.'); end

    config.steppedsine.freq.Start               = p(3);
    config.steppedsine.freq.End                 = p(4);
    config.steppedsine.freq.Step                = p(5);
    config.steppedsine.sampling.Rate            = p(6);
    config.steppedsine.sampling.SamplesPerBlock = p(7);
    config.steppedsine.system.BeamType          = p(8);
    config.steppedsine.sampling.Transient       = 0;
    config.steppedsine.rows.N                   = height(config.steppedsine.data.Table);

    config.steppedsine.NumOfBlocks = floor(config.steppedsine.rows.N / config.steppedsine.sampling.SamplesPerBlock);

    % Amplitude settings
    config.steppedsine.amplitude.Min    = 0.5;
    config.steppedsine.amplitude.Max    = 8;
    config.steppedsine.amplitude.Step   = 0.5;
    config.steppedsine.amplitude.Vector = 0.5:0.5:8;
    config.steppedsine.amplitude.N      = numel(config.steppedsine.amplitude.Vector);
    Na = config.steppedsine.amplitude.N;
    config.steppedsine.voltage.Vector = linspace(8/Na, 8, Na);

    if mod(config.steppedsine.rows.N, Na) ~= 0
        error('Total rows not divisible by number of amplitude blocks.');
    end
    config.steppedsine.rows.PerAmp = config.steppedsine.rows.N / Na;
    config.steppedsine.freq.N      = round(config.steppedsine.rows.PerAmp / config.steppedsine.sampling.SamplesPerBlock);

    % System-dependent settings
    switch config.steppedsine.system.BeamType
        case 1; config.steppedsine.system.Mass = 0.71345;
        case 2; config.steppedsine.system.Mass = 0.22821;
        case 3; config.steppedsine.system.Mass = 0.23304;
        otherwise; error('Incorrect beam type.');
    end
    config.steppedsine.channels.ResponseStart = 3;
    config.steppedsine.channels.ResponseEnd   = 6;
    config.steppedsine.channels.Excitation    = 2;
    config.steppedsine.channels.Response      = 3:6;

    const.g  = 9.80665;
    const.fc = 35;

    Nf = config.steppedsine.freq.N;

    % Per-temperature preallocation
    proc.FRF.X    = NaN(Na, Nf);
    proc.FRF.Y    = NaN(Na, Nf);
    proc.FRF.H    = NaN(Na, Nf);
    proc.FRF.aMag = NaN(Na, Nf);
    proc.FRF.dMag = NaN(Na, Nf);
    proc.FRF.HMag    = NaN(Na, Nf);
    proc.FRF.HMagS   = NaN(Na, Nf);
    proc.FRF.HPhase  = NaN(Na, Nf);
    proc.FRF.HPhaseS = NaN(Na, Nf);

    proc.peaks.ResponseAcc = NaN(Na, Nf);
    proc.energy.Dissipated = NaN(Na, Nf);

    proc.amp = repmat(struct( ...
        'Voltage',       [], ...
        'PeakAmp',       NaN, ...
        'Zeta',          NaN, ...
        'DampingMethod', '', ...
        'HPB',           struct('Low',NaN,'High',NaN), ...
        'Circle',        struct('centre',NaN,'radius',NaN,'quality',NaN,'nPoints',0), ...
        'Lean',          struct('Left',NaN,'Right',NaN,'Asym',NaN), ...
        'Energy',        struct('Mean',NaN,'EDR',NaN,'NDI',NaN), ...
        'Fit',           struct('fDense',[],'yFit',[]), ...
        'Stiffness',     struct('k_eff',NaN,'k_ref',NaN,'deltaK',NaN)), Na, 1);

    backbone.response.f = NaN(1, Na);
    backbone.response.y = NaN(1, Na);

    proc.stiffness.f_res  = NaN(Na, 1);
    proc.stiffness.k_eff  = NaN(Na, 1);
    proc.stiffness.deltaK = NaN(Na, 1);

    fTest = (config.steppedsine.freq.Start : config.steppedsine.freq.Step : ...
             config.steppedsine.freq.Start + (Nf-1)*config.steppedsine.freq.Step)';

    % Detect descending sweep, store direction, then normalise Start/End
    sweepDir = sign(mean(diff( ...
        config.steppedsine.data.Table{1:config.steppedsine.sampling.SamplesPerBlock:end,1})));
    if sweepDir < 0
        blockOrder = reshape(1:config.steppedsine.rows.N, config.steppedsine.sampling.SamplesPerBlock, []);
        flippedIdx = blockOrder(:, end:-1:1); flippedIdx = flippedIdx(:);
        config.steppedsine.data.Table = config.steppedsine.data.Table(flippedIdx, :);
        fTest = flipud(fTest);
        % Normalise so Start < End (fixes xlim and any other range-dependent code)
        freqLims = sort([config.steppedsine.freq.Start, config.steppedsine.freq.End]);
        config.steppedsine.freq.Start = freqLims(1);
        config.steppedsine.freq.End   = freqLims(2);
    end
    % Store sweep direction for hysteresis analysis (+1 = up, -1 = down)
    config.steppedsine.sweep.direction = sweepDir;

    % Retrieve impact-test reference values.
    % k_ref_impact comes from the system-specific impact test (correct — each
    % system's own baseline stiffness at T=20°C).
    % m_eff uses the MONOLITHIC physical mass for both systems so that k_eff,
    % c_eq and all regression coefficients are on the same absolute scale and
    % the jointed system correctly appears softer than the monolithic.
    if strcmp(config.system.label, 'Monolithic')
        k_ref_impact = proc.impact.monolithic.stiffness;
    else
        k_ref_impact = proc.impact.jointed.stiffness;
    end
    m_eff = proc.impact.monolithic.m_eff;   % shared physical mass — same for both systems
    T_impact  = 20;
    T_current = config.steppedsine.temperature.Vector(Z);

% =========================================================================
% AMPLITUDE LOOP
% =========================================================================

    for a = 1:Na

        proc.amp(a).Voltage = config.steppedsine.voltage.Vector(a);
        aStart = (a-1) * config.steppedsine.rows.PerAmp + 1;

        % Frequency loop — signal processing and FRF
        for k = 1:Nf
            r0 = aStart + (k-1)*config.steppedsine.sampling.SamplesPerBlock;
            r1 = r0 + config.steppedsine.sampling.SamplesPerBlock - 1;

            t   = config.steppedsine.data.Table{r0:r1, 1};
            t   = t - t(1);
            exc = detrend(config.steppedsine.data.Table{r0:r1, config.steppedsine.channels.Excitation});
            res = detrend(config.steppedsine.data.Table{r0:r1, config.steppedsine.channels.Response});
            res = mean(res, 2);

            fs = 1 / mean(diff(t));
            [bB, aB] = butter(2, const.fc/(fs/2), 'high');
            res = filtfilt(bB, aB, res);

            nt  = round(config.steppedsine.sampling.Transient * config.steppedsine.sampling.SamplesPerBlock);
            t   = t(nt+1:end);
            exc = exc(nt+1:end);
            res = res(nt+1:end);

            vel  = cumtrapz(t, res);
            vel  = vel - polyval(polyfit(t,vel,1), t);
            dis  = cumtrapz(t, vel);
            dis  = dis - polyval(polyfit(t,dis,2), t);

            proc.energy.Dissipated(a,k) = trapz(t, exc .* vel);

            fExc  = fTest(k);
            omega = 2*pi*fExc;

            proc.FRF.X(a,k) = 2/numel(exc) * sum(exc .* exp(-1i*2*pi*fExc*t));
            proc.FRF.Y(a,k) = 2/numel(res) * sum(res .* exp(-1i*2*pi*fExc*t));
            % Pipeline-normalised complex response quantity. Because the denominator
            % includes the specimen-mass and g scaling used in the acquisition workflow,
            % this is not automatically interpreted as a calibrated force/acceleration FRF.
            proc.FRF.H(a,k) = proc.FRF.Y(a,k) / (config.steppedsine.system.Mass * const.g * proc.FRF.X(a,k));

            proc.FRF.aMag(a,k) = abs(proc.FRF.Y(a,k));
            proc.FRF.dMag(a,k) = proc.FRF.aMag(a,k) / omega^2;

            proc.peaks.ResponseAcc(a,k) = localPeakEstimate(res);
        end

        proc.FRF.HMag(a,:)    = abs(proc.FRF.H(a,:));
        proc.FRF.HMagS(a,:)   = smoothdata(proc.FRF.HMag(a,:),  'movmean', 3);
        proc.FRF.HPhase(a,:)  = angle(proc.FRF.H(a,:));
        proc.FRF.HPhaseS(a,:) = smoothdata(proc.FRF.HPhase(a,:),'movmean', 3);

        % Resonance peak — PCHIP on acceleration response
        y = proc.peaks.ResponseAcc(a,:);
        proc.amp(a).Fit.fRaw = fTest(:);
        proc.amp(a).Fit.yRaw = y(:);
        if all(isnan(y)) || all(y == 0); continue; end

        validPts = ~isnan(y);
        pp = pchip(fTest(validPts), y(validPts));
        fD = linspace(min(fTest(validPts)), max(fTest(validPts)), ...
                      max(200, round((max(fTest(validPts))-min(fTest(validPts)))/0.001)));
        yF = ppval(pp, fD);

        [yPk, iPk] = max(yF);
        fPk        = fD(iPk);
        omega_pk   = 2*pi*fPk;

        proc.amp(a).PeakAmp   = yPk;
        peakAmpMatrix(a, Z)   = yPk;
        peakFreqMatrix(a, Z)  = fPk;
        proc.amp(a).Fit.fDense = fD;
        proc.amp(a).Fit.yFit   = yF;
        proc.amp(a).Energy.Mean = mean(abs(proc.energy.Dissipated(a,:)));

        % Effective stiffness
        if exist('c1_fit','var') && ~isnan(c1_fit)
            k_ref_T = k_ref_impact * (1 + c1_fit * (T_current - T_impact));
        else
            k_ref_T = k_ref_impact;
        end
        k_eff_loc = m_eff * omega_pk^2;

        proc.amp(a).Stiffness.k_eff  = k_eff_loc;
        proc.amp(a).Stiffness.k_ref  = k_ref_T;
        proc.amp(a).Stiffness.deltaK = (k_eff_loc - k_ref_T) / k_ref_T;
        proc.stiffness.f_res(a)  = fPk;
        proc.stiffness.k_eff(a)  = k_eff_loc;
        proc.stiffness.deltaK(a) = proc.amp(a).Stiffness.deltaK;

        % Harmonic displacement estimate reconstructed from acceleration at fPk
        HdMag = proc.FRF.dMag(a, :);
        validH_d = ~isnan(HdMag);
        if sum(validH_d) >= 3
            dMag_atPk = ppval(pchip(fTest(validH_d), HdMag(validH_d)), fPk);
        else
            [~, iPk_snap] = min(abs(fTest - fPk));
            dMag_atPk = HdMag(iPk_snap);
        end
        A_metres = dMag_atPk * const.g;
        fprintf('A_metres (a=%d, Z=%d): %.6e m  (%.4f um)\n', a, Z, A_metres, A_metres*1e6);
        AMatrix(a, Z)    = A_metres;
        kEffMatrix(a, Z) = k_eff_loc;

        % Damping — Nyquist circle fit, fallback to half-power bandwidth
        Hc = proc.FRF.H(a, :);
        HmagFull = abs(Hc);
        validHm  = ~isnan(HmagFull) & HmagFull > 0;
        if sum(validHm) >= 3
            HmFit  = ppval(pchip(fTest(validHm), HmagFull(validHm)), fD);
            aboveH = HmFit >= max(HmFit)/sqrt(2);
            if any(aboveH)
                bw_est = fD(find(aboveH,1,'last')) - fD(find(aboveH,1,'first'));
                bw_est = max(bw_est, 4*config.steppedsine.freq.Step);
            else
                bw_est = 4*config.steppedsine.freq.Step;
            end
        else
            bw_est = 4*config.steppedsine.freq.Step;
        end

        bw_factor = 2.5;
        win_mask  = fTest' >= (fPk - bw_factor*bw_est) & ...
                    fTest' <= (fPk + bw_factor*bw_est) & ...
                    ~isnan(Hc) & abs(Hc) > 0;
        Hc_win = Hc(win_mask);
        f_win  = fTest(win_mask);

        zeta_circle = NaN;
        circle_ok   = false;

        if numel(Hc_win) >= 5
            xp = real(Hc_win(:));
            yp = imag(Hc_win(:));
            coef = [xp, yp, ones(size(xp))] \ -(xp.^2 + yp.^2);
            xc = -coef(1)/2;   yc = -coef(2)/2;
            R2c = xc^2 + yc^2 - coef(3);
            if R2c > 0
                Rc = sqrt(R2c);
                dist_from_circle = abs(sqrt((xp-xc).^2 + (yp-yc).^2) - Rc);
                circle_quality   = mean(dist_from_circle) / (Rc + eps);

                [~, i_near] = min(abs(f_win - fPk));
                i_lo2 = max(1, i_near-1);
                i_hi2 = min(numel(f_win), i_near+1);

                if i_lo2 ~= i_hi2 && circle_quality < 0.3
                    dtheta = abs(angle(Hc_win(i_hi2) - complex(xc,yc)) - ...
                                 angle(Hc_win(i_lo2) - complex(xc,yc)));
                    if dtheta > pi; dtheta = 2*pi - dtheta; end
                    domega = 2*pi*(f_win(i_hi2) - f_win(i_lo2));
                    if domega > 0 && dtheta > 0
                        zeta_circle = 1 / (omega_pk * (dtheta/domega));
                        circle_ok = zeta_circle > 0.001 && zeta_circle < 0.15;
                    end
                end

                proc.amp(a).Circle.centre  = complex(xc, yc);
                proc.amp(a).Circle.radius  = Rc;
                proc.amp(a).Circle.quality = circle_quality;
                proc.amp(a).Circle.nPoints = numel(Hc_win);
            end
        end

        if circle_ok
            proc.amp(a).Zeta          = zeta_circle;
            proc.amp(a).DampingMethod = 'CircleFit';
        else
            yH = yPk / sqrt(2);
            f1 = findHalfPowerIntersection(fD, yF, iPk, 'left',  yH);
            f2 = findHalfPowerIntersection(fD, yF, iPk, 'right', yH);
            if ~isnan(f1) && ~isnan(f2) && f2 > f1
                proc.amp(a).Zeta          = (f2 - f1) / (2*fPk);
                proc.amp(a).HPB.Low       = f1;
                proc.amp(a).HPB.High      = f2;
                proc.amp(a).DampingMethod = 'HalfPower';
            else
                proc.amp(a).Zeta          = NaN;
                proc.amp(a).DampingMethod = 'Failed';
            end
        end

        % Peak asymmetry (always from PCHIP)
        yH = yPk / sqrt(2);
        f1l = findHalfPowerIntersection(fD, yF, iPk, 'left',  yH);
        f2l = findHalfPowerIntersection(fD, yF, iPk, 'right', yH);
        if ~isnan(f1l) && ~isnan(f2l) && f2l > f1l
            proc.amp(a).Lean.Left  = (fPk - f1l) / (f2l - f1l);
            proc.amp(a).Lean.Right = (f2l - fPk) / (f2l - f1l);
            proc.amp(a).Lean.Asym  = (f2l + f1l - 2*fPk) / (f2l - f1l);
        end

        backbone.response.f(a) = fPk;
        backbone.response.y(a) = yPk;

        zetaMatrix(a, Z) = proc.amp(a).Zeta;
        if ~isnan(proc.amp(a).Zeta)
            cEqMatrix(a, Z) = 2 * m_eff * proc.amp(a).Zeta * omega_pk;
        end

    end  % amplitude loop

    % =========================================================================
% EXTRACT EXCITATION ACCELERATION LEVELS (g RMS and peak)
% =========================================================================

exc_g_rms  = NaN(Na, 1);
exc_g_peak = NaN(Na, 1);

for a = 1:Na
    aStart = (a-1) * config.steppedsine.rows.PerAmp + 1;
    aEnd   = aStart + config.steppedsine.rows.PerAmp - 1;

    % Pull the full excitation block for this amplitude (all frequencies)
    exc_raw = config.steppedsine.data.Table{aStart:aEnd, ...
                  config.steppedsine.channels.Excitation};
    exc_raw = detrend(exc_raw);

    % Convert from raw units to g
    % (your data is in m/s^2 based on const.g usage — divide to get g)
    exc_g = exc_raw;% / const.g;

    exc_g_rms(a)  = rms(exc_g);
    exc_g_peak(a) = max(abs(exc_g));
end

fprintf('\nExcitation levels at T = %d°C:\n', T_current);
fprintf('  Amp level | RMS (g) | Peak (g)\n');
for a = 1:Na
    fprintf('  %9.1f | %7.4f | %8.4f\n', ...
        config.steppedsine.amplitude.Vector(a), exc_g_rms(a), exc_g_peak(a));
end

    fprintf('...done.\n\n');

    % Legacy uncalibrated energy-proxy normalisation (not a physical dissipation measurement)
    E0 = proc.amp(1).Energy.Mean;
    for a = 1:Na
        proc.amp(a).Energy.EDR = proc.amp(a).Energy.Mean / (E0 + eps);
        if ~isnan(proc.amp(a).PeakAmp) && proc.amp(a).PeakAmp > 0
            proc.amp(a).Energy.NDI = proc.amp(a).Energy.Mean / proc.amp(a).PeakAmp^2;
        end
        proc.amp(a).Voltage = config.steppedsine.voltage.Vector(a);
    end

    % Backbone polynomial fit (quadratic max)
    validBB = ~isnan(backbone.response.f) & ~isnan(backbone.response.y);
    if sum(validBB) >= 2
        aIdx   = find(validBB);
        pOrder = min(2, numel(aIdx)-1);
        p_fx   = polyfit(aIdx, backbone.response.f(aIdx), pOrder);
        p_fy   = polyfit(aIdx, backbone.response.y(aIdx), pOrder);
        aDense = linspace(1, max(aIdx), 10000);
        fCurve = polyval(p_fx, aDense);
        yCurve = polyval(p_fy, aDense);
        if sum(diff(sign(diff(fCurve))) ~= 0) > 1
            fCurve = polyval(polyfit(aIdx, backbone.response.f(aIdx), 1), aDense);
        end
    else
        fCurve = NaN;  yCurve = NaN;
    end

    proc.steppedsine.backbone{Z} = struct( ...
        'f',     fCurve(:), ...
        'y',     yCurve(:), ...
        'f_pts', backbone.response.f(:), ...
        'y_pts', backbone.response.y(:), ...
        'A_pts', AMatrix(:, Z), ...
        'k_pts', kEffMatrix(:, Z));
    proc.steppedsine.amp{Z} = proc.amp;
    proc.steppedsine.FRF{Z} = struct('HMagS', proc.FRF.HMagS, ...
                                      'HMag',  proc.FRF.HMag, ...
                                      'f',     fTest);

    % Legacy descriptive onset estimate
    dK     = proc.stiffness.deltaK;
    valid  = ~isnan(dK);
    onsetIdx = NaN;
    if sum(valid) >= 2
        exceed = abs(dK)*100 > onset_thresh_pct & valid;
        for a = 1:Na-1
            if exceed(a) && exceed(a+1)
                onsetIdx = a; break;
            end
        end
    end
    if isnan(onsetIdx)
        summaryOnset(Z)   = config.steppedsine.amplitude.Vector(Na);
        summaryOnsetUm(Z) = AMatrix(Na, Z) * 1e6;
    else
        summaryOnset(Z)   = config.steppedsine.amplitude.Vector(onsetIdx);
        summaryOnsetUm(Z) = AMatrix(onsetIdx, Z) * 1e6;
    end

    % Provisional descriptive label; replaced later by BFS classification
    dK_valid = dK(valid);
    summaryNonlinearityType{Z} = 'Softening';
    if mean(dK_valid) > 0; summaryNonlinearityType{Z} = 'Hardening'; end

end  % temperature loop

% =========================================================================
% DUFFING COEFFICIENT IDENTIFICATION
% f_res^2 = f0^2 + alpha * A^2  (cubic Duffing backbone linearisation)
% alpha > 0 = hardening;  alpha < 0 = softening
% =========================================================================

T_vec = config.steppedsine.temperature.Vector(:);
nT = config.steppedsine.nTemps;
duffing_alpha = NaN(1, nT);
duffing_f0sq  = NaN(1, nT);
duffing_R2    = NaN(1, nT);

for Z = 1:nT
    f_pts = proc.steppedsine.backbone{Z}.f_pts(:);   % resonance freq per amplitude
    A_pts = proc.steppedsine.backbone{Z}.A_pts(:);   % displacement per amplitude [m]
    valid = ~isnan(f_pts) & ~isnan(A_pts) & A_pts > 0;
    if sum(valid) < 4; continue; end

    f2 = f_pts(valid).^2;
    A2 = (A_pts(valid)*1e6).^2;   % convert to µm^2 for numerically stable fit

    X    = [ones(sum(valid),1), A2];
    beta = X \ f2;
    duffing_f0sq(Z)  = beta(1);
    duffing_alpha(Z) = beta(2);   % Hz^2 / µm^2
    fprintf('Duffing fit Z=%d: alpha=%.6f Hz2/um2  f0=%.4f Hz\n', Z, beta(2), sqrt(beta(1)));

    f2_pred = X * beta;
    ss_res  = sum((f2 - f2_pred).^2);
    ss_tot  = sum((f2 - mean(f2)).^2);
    if ss_tot > 0; duffing_R2(Z) = 1 - ss_res/ss_tot; end
end

proc.thermo.duffing.alpha  = duffing_alpha;
proc.thermo.duffing.f0sq   = duffing_f0sq;
proc.thermo.duffing.R2     = duffing_R2;
proc.thermo.duffing.T_vec  = T_vec(:)';

fprintf('Duffing alpha(T): min=%.3f  max=%.3f  Hz^2/um^2\n', ...
    min(duffing_alpha,[],'omitnan'), max(duffing_alpha,[],'omitnan'));
fprintf('R2 range: %.3f to %.3f\n', ...
    min(duffing_R2,[],'omitnan'), max(duffing_R2,[],'omitnan'));

% =========================================================================
% ALPHA(T) POLYNOMIAL FIT — Global Thermal-Mechanical Duffing Model
% f_res^2(A,T) = f0^2(T) + alpha(T) * A^2
% where alpha(T) = sum_k p_k * T^k  (polynomial in absolute temperature)
% and   f0^2(T)  = sum_k q_k * T^k
% =========================================================================

T_abs   = proc.thermo.duffing.T_vec(:);       % absolute temperature [°C]
alpha_v = proc.thermo.duffing.alpha(:);        % Hz^2/µm^2 per temperature
f0sq_v  = proc.thermo.duffing.f0sq(:);        % Hz^2 per temperature
R2_v    = proc.thermo.duffing.R2(:);

% --- Degree selection: fit poly orders 1–4, pick by adjusted R² ---
% Weight by per-temperature Duffing R² so unreliable alpha estimates
% (cold jointed, R²~0.03) contribute less to the polynomial fit.
w = max(R2_v, 0);   % zero-floor negative R² values

alpha_poly = struct();
f0sq_poly  = struct();

for deg = 1:4
    % Alpha(T) polynomial
    % Only positive-weight, finite observations contribute to weighted fitting.
    valid_a = isfinite(alpha_v) & isfinite(T_abs) & isfinite(w) & w > 0;
    n_alpha = nnz(valid_a);
    n_parameters = deg + 1;  % degree d polynomial includes the intercept
    if n_alpha <= n_parameters
        error('research_pipeline:InsufficientDuffingData', ...
            'Degree %d alpha(T) model needs >%d positive-weight observations; found %d.', ...
            deg, n_parameters, n_alpha);
    end
    p_a = wpolyfit(T_abs(valid_a), alpha_v(valid_a), w(valid_a), deg);
    alpha_pred = polyval(p_a, T_abs(valid_a));
    w_alpha = w(valid_a);
    y_alpha = alpha_v(valid_a);
    mean_alpha = sum(w_alpha .* y_alpha) / sum(w_alpha);
    ss_res = sum(w_alpha .* (y_alpha-alpha_pred).^2);
    ss_tot = sum(w_alpha .* (y_alpha-mean_alpha).^2);
    R2_a = 1 - ss_res / max(ss_tot, eps);
    % Adjusted R2 denominator is n - p, where p includes the intercept.
    adjR2_a = 1 - (1-R2_a)*(n_alpha-1)/(n_alpha-n_parameters);

    % f0sq(T) uses its own valid sample count and unweighted fit.
    valid_f = isfinite(f0sq_v) & isfinite(T_abs);
    n_f0 = nnz(valid_f);
    if n_f0 <= n_parameters
        error('research_pipeline:InsufficientDuffingData', ...
            'Degree %d f0sq(T) model needs >%d observations; found %d.', ...
            deg, n_parameters, n_f0);
    end
    p_f = polyfit(T_abs(valid_f), f0sq_v(valid_f), deg);
    f0sq_pred = polyval(p_f, T_abs(valid_f));
    ss_res_f = sum((f0sq_v(valid_f)-f0sq_pred).^2);
    ss_tot_f = sum((f0sq_v(valid_f)-mean(f0sq_v(valid_f))).^2);
    R2_f = 1 - ss_res_f / max(ss_tot_f, eps);
    adjR2_f = 1 - (1-R2_f)*(n_f0-1)/(n_f0-n_parameters);

    alpha_poly(deg).coeffs = p_a;
    alpha_poly(deg).R2     = R2_a;
    alpha_poly(deg).adjR2  = adjR2_a;
    f0sq_poly(deg).coeffs  = p_f;
    f0sq_poly(deg).R2      = R2_f;
    f0sq_poly(deg).adjR2   = adjR2_f;

    fprintf('deg %d | alpha(T): R²=%.4f adjR²=%.4f | f0²(T): R²=%.4f adjR²=%.4f\n', ...
        deg, R2_a, adjR2_a, R2_f, adjR2_f);
end

% --- Select best degree by adjusted R² for alpha(T) ---
adjR2_vals = [alpha_poly.adjR2];
[~, best_deg] = max(adjR2_vals);
p_alpha_best  = alpha_poly(best_deg).coeffs;
p_f0sq_best   = f0sq_poly(best_deg).coeffs;
fprintf('\nSelected polynomial degree: %d  (alpha adjR²=%.4f)\n', ...
    best_deg, alpha_poly(best_deg).adjR2);

% --- Global surface prediction: f_res^2(A,T) ---
% Predict over a dense (T, A) grid for plotting and export
T_dense = linspace(min(T_abs), max(T_abs), 200);
A_dense = linspace(0, max(AMatrix(:), [], 'omitnan')*1e6, 100);   % µm
[T_grid, A_grid] = meshgrid(T_dense, A_dense);

alpha_grid = polyval(p_alpha_best, T_grid);
f0sq_grid  = polyval(p_f0sq_best,  T_grid);
fres_grid  = sqrt(max(f0sq_grid + alpha_grid .* A_grid.^2, 0));

% --- Store results ---
proc.thermo.duffing.alpha_poly = alpha_poly;
proc.thermo.duffing.f0sq_poly  = f0sq_poly;
proc.thermo.duffing.best_deg   = best_deg;
proc.thermo.duffing.p_alpha    = p_alpha_best;
proc.thermo.duffing.p_f0sq     = p_f0sq_best;
proc.thermo.duffing.T_grid     = T_grid;
proc.thermo.duffing.A_grid     = A_grid;
proc.thermo.duffing.fres_grid  = fres_grid;

% --- Diagnostic figure ---
figure('Name', sprintf('Alpha(T) Polynomial — %s', config.system.label));

subplot(1,3,1);
scatter(T_abs(valid_a), alpha_v(valid_a), 60, w(valid_a), 'filled'); hold on;
T_plot = linspace(min(T_abs), max(T_abs), 300);
plot(T_plot, polyval(p_alpha_best, T_plot), 'r-', 'LineWidth', 2);
yline(0, 'k--');
xlabel('Temperature (°C)'); ylabel('\alpha (Hz^2/\mum^2)');
title(sprintf('\\alpha(T) — degree %d  adjR^2=%.3f', best_deg, alpha_poly(best_deg).adjR2));
colorbar; clabel_text = 'Duffing R^2'; % colour = per-T reliability
grid on;

subplot(1,3,2);
scatter(T_abs(valid_f), f0sq_v(valid_f), 60, 'k', 'filled'); hold on;
plot(T_plot, polyval(p_f0sq_best, T_plot), 'r-', 'LineWidth', 2);
xlabel('Temperature (°C)'); ylabel('f_0^2 (Hz^2)');
title(sprintf('f_0^2(T) — degree %d  R^2=%.3f', best_deg, f0sq_poly(best_deg).R2));
grid on;

subplot(1,3,3);
surf(T_grid, A_grid, fres_grid, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
xlabel('Temperature (°C)'); ylabel('Amplitude (\mum)'); zlabel('f_{res} (Hz)');
title('Global Model: f_{res}(A,T)');
colorbar; view(45, 30); grid on;

% =========================================================================
% MINDLIN POWER-LAW FIT: delta_K ~ C * A^n  above the zero-crossing
% Uses Duffing f0 as reference, consistent with Table 2 deltaK definition
% =========================================================================

powerlaw_n   = NaN(1, nT);
powerlaw_C   = NaN(1, nT);
powerlaw_R2  = NaN(1, nT);

for Z = 1:nT
    f_pts = proc.steppedsine.backbone{Z}.f_pts(:);
    A_pts = proc.steppedsine.backbone{Z}.A_pts(:) * 1e6;  % µm

    valid = ~isnan(f_pts) & ~isnan(A_pts) & A_pts > 0;
    if sum(valid) < 4; continue; end

    % Use Duffing-extrapolated f0 as reference — same as Table 2 deltaK
    f0_sq = duffing_f0sq(Z);
    if isnan(f0_sq) || f0_sq <= 0; continue; end
    f0_T  = sqrt(f0_sq);

    dK_vec = ((f_pts(valid).^2 - f0_sq) / f0_sq) * 100;  % % deviation
    A_vec  = A_pts(valid);
    fprintf('dK Z=%d: A range %.4f-%.4f um  dK range %.3f-%.3f%%\n', Z, min(A_vec), max(A_vec), min(dK_vec), max(dK_vec));

    % Find crossover: dK transitions from negative to positive
    cross_idx = NaN;
    for ai = 1:numel(dK_vec)-1
        if dK_vec(ai) < 0 && dK_vec(ai+1) >= 0
            cross_idx = ai + 1; break;
        end
    end
    if isnan(cross_idx); continue; end

    A_above  = A_vec(cross_idx:end);
    dK_above = dK_vec(cross_idx:end);
    pos      = dK_above > 0;
    if sum(pos) < 3; continue; end

    logA  = log(A_above(pos));
    logdK = log(dK_above(pos));
    p_fit = polyfit(logA, logdK, 1);
    n_fit = p_fit(1);
    C_fit = exp(p_fit(2));

    dK_pred = C_fit * A_above(pos).^n_fit;
    ss_res  = sum((dK_above(pos) - dK_pred).^2);
    ss_tot  = sum((dK_above(pos) - mean(dK_above(pos))).^2);

    powerlaw_n(Z)  = n_fit;
    powerlaw_C(Z)  = C_fit;
    powerlaw_R2(Z) = 1 - ss_res / (ss_tot + eps);
end

proc.thermo.powerlaw.n   = powerlaw_n;
proc.thermo.powerlaw.C   = powerlaw_C;
proc.thermo.powerlaw.R2  = powerlaw_R2;

% Print table
fprintf('\nMindlin power-law exponent n(T):\n');
fprintf('%-10s %-8s %-8s %-8s\n', 'T (degC)', 'n', 'C', 'R2');
for Z = 1:nT
    T_val = config.steppedsine.temperature.Vector(Z);
    if ~isnan(powerlaw_n(Z))
        fprintf('%-10d %-8.3f %-8.4f %-8.4f\n', ...
            T_val, powerlaw_n(Z), powerlaw_C(Z), powerlaw_R2(Z));
    else
        fprintf('%-10d  — (no crossover)\n', T_val);
    end
end
fprintf('Mean n = %.3f,  std = %.3f\n', ...
    mean(powerlaw_n,'omitnan'), std(powerlaw_n,'omitnan'));

% =========================================================================
% REGRESSION
% =========================================================================

T_vec    = config.steppedsine.temperature.Vector(:);
T0       = 20;
dT_vec   = T_vec - T0;
nT       = numel(T_vec);

A_flat    = AMatrix(:);
dT_flat   = repmat(dT_vec', Na, 1);  dT_flat  = dT_flat(:);
T_flat    = repmat(T_vec',  Na, 1);  T_flat   = T_flat(:);
kEff_flat = kEffMatrix(:);
cEq_flat  = cEqMatrix(:);

valid_k     = ~isnan(A_flat) & ~isnan(kEff_flat);
valid_c_raw = ~isnan(A_flat) & ~isnan(cEq_flat);

A_v  = A_flat(valid_k);   dT_v = dT_flat(valid_k);   k_v = kEff_flat(valid_k);
A_c_all  = A_flat(valid_c_raw);
dT_c_all = dT_flat(valid_c_raw);
T_c_all  = T_flat(valid_c_raw);
c_v_all  = cEq_flat(valid_c_raw);

% -------------------------------------------------------------------------
% STIFFNESS-RELATED REGRESSION — linear-T piecewise
% Model: k = k_ref * (1 + c1*dT + [g0_regime + g1_regime*dT]*A^2)
% c2*dT^2 removed: degraded R² in both systems.
%
% T* selection: use fixed value when config.steppedsine.T_star_fixed is provided
% is set; otherwise sweep 0:10:60 and pick best R².
% -------------------------------------------------------------------------

use_fixed_Tstar = isfield(config.steppedsine,'T_star_fixed') && ...
                  ~isnan(config.steppedsine.T_star_fixed);

T_star_candidates = 0:10:60;
R2_k_sweep        = NaN(size(T_star_candidates));

if use_fixed_Tstar
    T_star_opt = config.steppedsine.T_star_fixed - T0;
    fprintf('\n=== T* fixed at %.0f°C (dT=%.0f) — R² sweep skipped ===\n', ...
        config.steppedsine.T_star_fixed, T_star_opt);
    r_fixed = double(dT_v >= T_star_opt);
    if sum(r_fixed<1)>=10 && sum(r_fixed>0)>=10
        Phi_f = [ones(size(A_v)), dT_v, ...
                 (1-r_fixed).*A_v.^2, (1-r_fixed).*dT_v.*A_v.^2, ...
                 r_fixed.*A_v.^2,      r_fixed.*dT_v.*A_v.^2];
        k_hat_f = Phi_f * (Phi_f \ k_v);
        R2_k = 1 - sum((k_v-k_hat_f).^2)/sum((k_v-mean(k_v)).^2);
    else
        R2_k = NaN;
        warning('Insufficient data in one regime for fixed T*=%.0f°C.', ...
            config.steppedsine.T_star_fixed);
    end
else
    for iStar = 1:numel(T_star_candidates)
        T_s = T_star_candidates(iStar);
        r_s = double(dT_v >= T_s);
        if sum(r_s<1)<10 || sum(r_s>0)<10; continue; end
        Phi_s = [ones(size(A_v)), dT_v, ...
                 (1-r_s).*A_v.^2, (1-r_s).*dT_v.*A_v.^2, ...
                 r_s.*A_v.^2,      r_s.*dT_v.*A_v.^2];
        k_hat = Phi_s * (Phi_s \ k_v);
        R2_k_sweep(iStar) = 1 - sum((k_v-k_hat).^2)/sum((k_v-mean(k_v)).^2);
    end
    [R2_k, iStar_best] = max(R2_k_sweep);
    T_star_opt = T_star_candidates(iStar_best);
end

regime_k = double(dT_v >= T_star_opt);

% Final linear-T piecewise fit
Phi_k  = [ones(size(A_v)), dT_v, ...
          (1-regime_k).*A_v.^2, (1-regime_k).*dT_v.*A_v.^2, ...
          regime_k.*A_v.^2,      regime_k.*dT_v.*A_v.^2];
theta_k = Phi_k \ k_v;
k_pred  = Phi_k * theta_k;

k_ref_fit   = theta_k(1);
c1_fit      = theta_k(2) / k_ref_fit;
gamma0_cold = theta_k(3) / k_ref_fit;
gamma1_cold = theta_k(4) / k_ref_fit;
gamma0_warm = theta_k(5) / k_ref_fit;
gamma1_warm = theta_k(6) / k_ref_fit;

% Global model (no piecewise split)
Phi_k_g    = [ones(size(A_v)), dT_v, A_v.^2, dT_v.*A_v.^2];
theta_k_g  = Phi_k_g \ k_v;
R2_k_global = 1 - sum((k_v - Phi_k_g*theta_k_g).^2) / sum((k_v-mean(k_v)).^2);
gamma0_fit  = theta_k_g(3)/theta_k_g(1);
gamma1_fit  = theta_k_g(4)/theta_k_g(1);

% Stiffness-related residual structure by temperature and amplitude index
k_resid = k_v - k_pred;
T_resid_mean = NaN(1,nT);  T_resid_std = NaN(1,nT);
A_resid_mean = NaN(1,Na);  A_resid_std = NaN(1,Na);
for iT_r = 1:nT
    idx_r = dT_v == dT_vec(iT_r);
    if sum(idx_r)>=2
        T_resid_mean(iT_r) = mean(k_resid(idx_r));
        T_resid_std(iT_r)  = std(k_resid(idx_r));
    end
end
A_vec_flat_idx = repmat((1:Na)',nT,1);
A_v_idx        = A_vec_flat_idx(valid_k);
for iA_r = 1:Na
    idx_r = A_v_idx == iA_r;
    if sum(idx_r)>=2
        A_resid_mean(iA_r) = mean(k_resid(idx_r));
        A_resid_std(iA_r)  = std(k_resid(idx_r));
    end
end

if strcmp(config.system.label,'Monolithic')
    f_ref_Hz = fM0;
else
    f_ref_Hz = fJ0;
end

fprintf('\n=== STIFFNESS (piecewise, linear-T) ===\n');
if use_fixed_Tstar
    fprintf('T_star=%.0f°C (FIXED)  k_ref=%.2f N/m  f_ref=%.4f Hz\n', T0+T_star_opt, k_ref_fit, f_ref_Hz);
else
    fprintf('T_star=%.0f°C (swept)  k_ref=%.2f N/m  f_ref=%.4f Hz\n', T0+T_star_opt, k_ref_fit, f_ref_Hz);
end
fprintf('  c1=%.4e /°C\n', c1_fit);
f_ratio_sq = (fJ0/fM0)^2;
fprintf('Sanity: (f_J/f_M)^2=%.4f  joint=%.1f%% of mono\n', f_ratio_sq, 100*f_ratio_sq);
fprintf('Cold: g0=%.4e  g1=%.4e  /m²\n', gamma0_cold, gamma1_cold);
fprintf('Warm: g0=%.4e  g1=%.4e  /m²\n', gamma0_warm, gamma1_warm);
fprintf('R²: piecewise=%.4f  global=%.4f\n', R2_k, R2_k_global);
[~,iT_maxres] = max(abs(T_resid_mean),[],'omitnan');
fprintf('Max |residual mean| by T: %.1f N/m at T=%.0f°C\n', ...
    abs(T_resid_mean(iT_maxres)), T_vec(iT_maxres));

% -------------------------------------------------------------------------
% EQUIVALENT-DAMPING REGRESSION — bounded exponential model
% Model: c_eq = c_ref*exp(d1*dT) + d2*A + d3*A^2
% Two-stage: (1) log-scale regression for d1 from low-amplitude data;
% (2) amplitude residuals for d2, d3.
% Hard bound: |d1| <= d1_max = 0.05 /°C.
% -------------------------------------------------------------------------

MAD_thresh     = 3.0;
c_outlier_flag = false(size(c_v_all));
for iT = 1:numel(unique(T_c_all))
    T_lvl = unique(T_c_all); T_lvl = T_lvl(iT);
    idx_T = T_c_all == T_lvl;
    c_T   = c_v_all(idx_T);
    c_med = median(c_T);
    mad_T = median(abs(c_T - c_med));
    if mad_T > eps
        c_outlier_flag(idx_T) = abs(c_T - c_med) > MAD_thresh * mad_T;
    end
end
n_flagged = sum(c_outlier_flag);
fprintf('\n=== DAMPING (MAD outliers: %d removed of %d) ===\n', n_flagged, numel(c_v_all));

valid_c = ~c_outlier_flag;
A_c  = A_c_all(valid_c);
dT_c = dT_c_all(valid_c);
T_c  = T_c_all(valid_c);
c_v  = c_v_all(valid_c);

isMonolithic = strcmp(config.system.label, 'Monolithic');

% Baseline comparison (unchanged model, for R² tracking)
Phi_c_old  = [ones(size(A_c)), dT_c, A_c, dT_c.*A_c];
c_pred_old = Phi_c_old * (Phi_c_old \ c_v);
R2_c_old   = 1 - sum((c_v-c_pred_old).^2)/sum((c_v-mean(c_v)).^2);

regime_c = double(T_c >= T_star_opt + T0);
cold_c   = 1 - regime_c;
warm_c   =     regime_c;

d1_max = 0.05;   % hard bound on |d1| in /°C


if isMonolithic
    [log_c_ref_cold, d1_cold] = fitLogDamping(dT_c, c_v, d1_max);
    c_ref_cold = exp(log_c_ref_cold);

    c_temp     = c_ref_cold * exp(d1_cold * dT_c);
    Phi_amp    = [A_c, A_c.^2];
    b_amp      = Phi_amp \ (c_v - c_temp);
    d2_cold = b_amp(1);   d3_cold = b_amp(2);
    c_ref_warm = c_ref_cold;  d1_warm = d1_cold;
    d2_warm    = d2_cold;     d3_warm = d3_cold;

    c_pred_pw = c_ref_cold*exp(d1_cold*dT_c) + d2_cold*A_c + d3_cold*A_c.^2;
    R2_c_pw   = 1 - sum((c_v-c_pred_pw).^2)/sum((c_v-mean(c_v)).^2);
    R2_c_v3   = R2_c_pw;

    fprintf('Single-regime: c_ref=%.4f  d1=%.4e /°C  d2=%.4e  d3=%.4e\n', ...
        c_ref_cold, d1_cold, d2_cold, d3_cold);
    fprintf('R²: old=%.4f  new(exp)=%.4f\n', R2_c_old, R2_c_pw);

else  % Jointed — piecewise regimes
    cold_idx = cold_c == 1;
    warm_idx = warm_c == 1;

    [log_c_ref_cold, d1_cold] = fitLogDamping(dT_c(cold_idx), c_v(cold_idx), d1_max);
    [log_c_ref_warm, d1_warm] = fitLogDamping(dT_c(warm_idx), c_v(warm_idx), d1_max);
    c_ref_cold = exp(log_c_ref_cold);
    c_ref_warm = exp(log_c_ref_warm);

    b_c = [A_c(cold_idx), A_c(cold_idx).^2] \ (c_v(cold_idx) - c_ref_cold*exp(d1_cold*dT_c(cold_idx)));
    b_w = [A_c(warm_idx), A_c(warm_idx).^2] \ (c_v(warm_idx) - c_ref_warm*exp(d1_warm*dT_c(warm_idx)));
    d2_cold = b_c(1);  d3_cold = b_c(2);
    d2_warm = b_w(1);  d3_warm = b_w(2);

    c_temp = zeros(size(c_v));
    c_temp(cold_idx) = c_ref_cold*exp(d1_cold*dT_c(cold_idx));
    c_temp(warm_idx) = c_ref_warm*exp(d1_warm*dT_c(warm_idx));
    c_pred_pw = c_temp + cold_c.*(d2_cold*A_c + d3_cold*A_c.^2) ...
                       + warm_c.*(d2_warm*A_c + d3_warm*A_c.^2);
    R2_c_pw = 1 - sum((c_v-c_pred_pw).^2)/sum((c_v-mean(c_v)).^2);

    Phi_c_v3  = [cold_c, cold_c.*dT_c, cold_c.*A_c, warm_c, warm_c.*dT_c, warm_c.*A_c];
    c_pred_v3 = Phi_c_v3 * (Phi_c_v3 \ c_v);
    R2_c_v3   = 1 - sum((c_v-c_pred_v3).^2)/sum((c_v-mean(c_v)).^2);

    fprintf('Cold: c_ref=%.4f  d1=%.4e /°C  d2=%.4e  d3=%.4e\n', c_ref_cold,d1_cold,d2_cold,d3_cold);
    fprintf('Warm: c_ref=%.4f  d1=%.4e /°C  d2=%.4e  d3=%.4e\n', c_ref_warm,d1_warm,d2_warm,d3_warm);
    fprintf('R²: old=%.4f  v3=%.4f  new(exp)=%.4f\n', R2_c_old, R2_c_v3, R2_c_pw);
end

d1_cold_n = d1_cold;   % /°C on log scale — no further normalisation
d2_cold_n = d2_cold/(c_ref_cold+eps);
d3_cold_n = d3_cold/(c_ref_cold+eps);
d1_warm_n = d1_warm;
d2_warm_n = d2_warm/(c_ref_warm+eps);
d3_warm_n = d3_warm/(c_ref_warm+eps);

% -------------------------------------------------------------------------
% BFS-BASED DESCRIPTIVE RESPONSE CLASSIFICATION
% Uses BFS sign (positive = hardening, negative = softening) from measured
% resonance frequencies. Threshold ±0.10% separates signal from noise.
% A secondary SNR check on deltaK retains Weak/Indeterminate for genuinely
% noisy temperature levels.
% -------------------------------------------------------------------------
BFS_thresh_pct = 0.10;

for Z = 1:nT
    f_pts   = proc.steppedsine.backbone{Z}.f_pts;
    valid_f = ~isnan(f_pts);

    if sum(valid_f) < 2
        summaryNonlinearityType{Z} = 'Indeterminate';
        continue;
    end

    f_valid = f_pts(valid_f);
    bfs_pct = (f_valid(end) - f_valid(1)) / f_valid(1) * 100;

    dK_vec   = arrayfun(@(a) proc.steppedsine.amp{Z}(a).Stiffness.deltaK, 1:Na);
    dK_clean = dK_vec(~isnan(dK_vec));
    dK_snr   = 0;
    if numel(dK_clean) >= 3
        dK_snr = abs(mean(dK_clean)) / (std(dK_clean) + eps);
    end

    if abs(bfs_pct) < BFS_thresh_pct && dK_snr < 3.0
        summaryNonlinearityType{Z} = 'Weak/Indeterminate';
    elseif bfs_pct >= BFS_thresh_pct
        summaryNonlinearityType{Z} = 'Hardening';
    elseif bfs_pct <= -BFS_thresh_pct
        summaryNonlinearityType{Z} = 'Softening';
    else
        summaryNonlinearityType{Z} = 'Weak/Indeterminate';
    end
end

% Post-hoc model-based onset
summaryOnset_model   = NaN(1, nT);
summaryOnsetUm_model = NaN(1, nT);
for Z = 1:nT
    g0 = gamma0_cold;
    if T_vec(Z)-T0 >= T_star_opt; g0 = gamma0_warm; end
    prev_exceed = false;
    for a = 1:Na
        A_m = AMatrix(a,Z);
        if ~isnan(A_m) && abs(g0 * A_m^2) * 100 > onset_thresh_pct
            if prev_exceed
                summaryOnset_model(Z)   = config.steppedsine.amplitude.Vector(a-1);
                summaryOnsetUm_model(Z) = AMatrix(a-1, Z) * 1e6;
                break;
            end
            prev_exceed = true;
        else
            prev_exceed = false;
        end
    end
    if isnan(summaryOnset_model(Z))
        summaryOnset_model(Z)   = config.steppedsine.amplitude.Vector(end);
        summaryOnsetUm_model(Z) = AMatrix(Na, Z) * 1e6;
    end
end

% Store results
proc.thermo.stiffness.f_ref_Hz = f_ref_Hz;
proc.thermo.stiffness.coeffs = struct('k_ref',k_ref_fit,'c1',c1_fit, ...
    'gamma0',gamma0_fit,'gamma1',gamma1_fit, ...
    'gamma0_cold',gamma0_cold,'gamma1_cold',gamma1_cold, ...
    'gamma0_warm',gamma0_warm,'gamma1_warm',gamma1_warm, ...
    'T_star_dT',T_star_opt,'T_star_abs',T0+T_star_opt, ...
    'T_star_fixed',use_fixed_Tstar);
proc.thermo.stiffness.R2              = R2_k;
proc.thermo.stiffness.R2_global       = R2_k_global;
proc.thermo.stiffness.R2_sweep        = R2_k_sweep;
proc.thermo.stiffness.T_star_sweep    = T_star_candidates;
proc.thermo.stiffness.A_data          = A_v;
proc.thermo.stiffness.dT_data         = dT_v;
proc.thermo.stiffness.k_meas          = k_v;
proc.thermo.stiffness.k_pred          = k_pred;
proc.thermo.stiffness.resid_by_T_mean = T_resid_mean;
proc.thermo.stiffness.resid_by_T_std  = T_resid_std;
proc.thermo.stiffness.resid_by_A_mean = A_resid_mean;
proc.thermo.stiffness.resid_by_A_std  = A_resid_std;

proc.thermo.damping.coeffs = struct('c_ref_cold',c_ref_cold,'c_ref_warm',c_ref_warm, ...
    'd1_cold',d1_cold,'d1_warm',d1_warm,'d2_cold',d2_cold,'d2_warm',d2_warm, ...
    'd3_cold',d3_cold,'d3_warm',d3_warm,'d1_cold_n',d1_cold_n,'d1_warm_n',d1_warm_n, ...
    'd2_cold_n',d2_cold_n,'d2_warm_n',d2_warm_n,'d3_cold_n',d3_cold_n,'d3_warm_n',d3_warm_n, ...
    'T_star_abs',T0+T_star_opt,'MAD_thresh',MAD_thresh,'n_outliers',n_flagged, ...
    'd1_bound',d1_max,'damping_model','exp_bounded');
proc.thermo.damping.R2           = R2_c_pw;
proc.thermo.damping.R2_v3        = R2_c_v3;
proc.thermo.damping.R2_old       = R2_c_old;
proc.thermo.damping.c_meas       = c_v;
proc.thermo.damping.c_pred       = c_pred_pw;
proc.thermo.damping.regime       = regime_c;
proc.thermo.damping.outlier_flag = c_outlier_flag;
proc.thermo.damping.c_meas_all   = c_v_all;
proc.thermo.damping.A_all        = A_c_all;
proc.thermo.damping.dT_all       = dT_c_all;
proc.thermo.onset.loop           = summaryOnset;
proc.thermo.onset.loop_um        = summaryOnsetUm;
proc.thermo.onset.model          = summaryOnset_model;
proc.thermo.onset.model_um       = summaryOnsetUm_model;

% =========================================================================
% CSV EXPORT
% =========================================================================

switch config.steppedsine.system.BeamType
    case 2;    systemLabel = 'Monolithic';
    case 3;    systemLabel = 'Jointed';
    otherwise; systemLabel = sprintf('BeamType%d', config.steppedsine.system.BeamType);
end

outlier_lookup = false(Na, nT);
for ai_ol = 1:Na
    for zi_ol = 1:nT
        flat_ol = (zi_ol-1)*Na + ai_ol;
        if flat_ol <= numel(valid_c_raw) && valid_c_raw(flat_ol)
            outlier_lookup(ai_ol,zi_ol) = c_outlier_flag(sum(valid_c_raw(1:flat_ol)));
        end
    end
end

exportRows = {};
for zi = 1:nT
    T_val  = config.steppedsine.temperature.Vector(zi);
    dT_val = T_val - T_vec(1);
    for ai = 1:Na
        exportRows{end+1} = {systemLabel, T_val, dT_val, ...
            config.steppedsine.amplitude.Vector(ai), AMatrix(ai,zi), kEffMatrix(ai,zi), ...
            (kEffMatrix(ai,zi)-k_ref_fit)/k_ref_fit, zetaMatrix(ai,zi), cEqMatrix(ai,zi), ...
            peakFreqMatrix(ai,zi), peakAmpMatrix(ai,zi), ...
            double(outlier_lookup(ai,zi)), proc.steppedsine.amp{zi}(ai).DampingMethod};
    end
end

obsTable = cell2table(vertcat(exportRows{:}), 'VariableNames', ...
    {'SystemType','Temperature_degC','DeltaT_degC','ExcitationAmplitude', ...
     'A_metres','k_eff_Npm','deltaK_normalised','Zeta','c_eq_Nspm', ...
     'ResonanceFreq_Hz','PeakAcc_g','DampingOutlierFlag','DampingMethod'});

% Export fields for the retained piecewise thermal model
coeffTable = table({systemLabel}, k_ref_fit, f_ref_Hz, c1_fit, ...
    R2_k, R2_k_global, T0+T_star_opt, double(use_fixed_Tstar), ...
    gamma0_cold, gamma1_cold, gamma0_warm, gamma1_warm, ...
    c_ref_cold, c_ref_warm, d1_cold, d1_warm, d2_cold, d2_warm, d3_cold, d3_warm, ...
    n_flagged, R2_c_pw, R2_c_v3, R2_c_old, ...
    'VariableNames', {'SystemType','k_ref_Npm','f_ref_Hz','c1_per_degC', ...
    'R2_stiffness_piecewise','R2_stiffness_global','T_star_degC','T_star_fixed', ...
    'gamma0_cold_per_m2','gamma1_cold_per_m2_degC','gamma0_warm_per_m2','gamma1_warm_per_m2_degC', ...
    'c_ref_cold_Nspm','c_ref_warm_Nspm','d1_cold_per_degC','d1_warm_per_degC', ...
    'd2_cold_Nspm_per_m','d2_warm_Nspm_per_m','d3_cold_Nspm_per_m2','d3_warm_Nspm_per_m2', ...
    'n_damping_outliers','R2_damping','R2_damping_v3','R2_damping_old'});

summaryRows = {};
for zi = 1:nT
    T_val = config.steppedsine.temperature.Vector(zi);
    bb    = proc.steppedsine.backbone{zi};
    fShift = NaN;
    if ~all(isnan(bb.f_pts))
        fpts = bb.f_pts(~isnan(bb.f_pts));
        if numel(fpts) >= 2; fShift = (fpts(end)-fpts(1))/fpts(1); end
    end
    zetas = arrayfun(@(a) proc.steppedsine.amp{zi}(a).Zeta, 1:Na);
    summaryRows{end+1} = {systemLabel, T_val, summaryOnset(zi), summaryOnsetUm(zi), ...
                           summaryNonlinearityType{zi}, fShift, mean(zetas,'omitnan')};
end
summaryTable = cell2table(vertcat(summaryRows{:}), 'VariableNames', ...
    {'SystemType','Temperature_degC','OnsetAmplitude_V','OnsetDisplacement_um', ...
     'NonlinearityType','BackboneFreqShift','MeanZeta'});

if ~exist('exportCache','var'); exportCache = struct(); end
exportCache.(systemLabel).obs     = obsTable;
exportCache.(systemLabel).coeffs  = coeffTable;
exportCache.(systemLabel).summary = summaryTable;

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
writetable(obsTable,     sprintf('export_%s_observations_%s.csv',  systemLabel, timestamp));
writetable(coeffTable,   sprintf('export_%s_coefficients_%s.csv',  systemLabel, timestamp));
writetable(summaryTable, sprintf('export_%s_summary_%s.csv',       systemLabel, timestamp));


% Power law export (jointed only — NaN rows included for completeness)
T_export = config.steppedsine.temperature.Vector(:);
powerlaw_table = table(T_export, powerlaw_n(:), powerlaw_C(:), powerlaw_R2(:), ...
    'VariableNames', {'Temperature_degC', 'Mindlin_n', 'Mindlin_C', 'Mindlin_R2'});
writetable(powerlaw_table, ...
    sprintf('export_%s_powerlaw_%s.csv', systemLabel, timestamp));
fprintf('Exported: power law n values — %d rows\n', height(powerlaw_table));

fprintf('Exported: %s obs=%d rows\n', systemLabel, height(obsTable));




% % =========================================================================
% % INTERACTIVE TEMPERATURE BROWSER
% % =========================================================================
%
% switch config.steppedsine.system.BeamType
%     case 2; figName = 'Temperature Browser — Monolithic';
%     case 3; figName = 'Temperature Browser — Jointed';
%     otherwise; figName = 'Temperature Browser';
% end
% fig = figure('Name', figName);
% ax  = axes(fig); hold(ax,'on'); grid(ax,'on');
% xlabel(ax,'Frequency (Hz)'); ylabel(ax,'Response Amplitude');
% set(ax,'YScale','log');
% xlim(sort([config.steppedsine.freq.Start, config.steppedsine.freq.End]));
% ylim([0 100]);
%
% setappdata(fig,'activeZ', 1);
% setappdata(fig,'temps',   config.steppedsine.temperature.Vector);
% setappdata(fig,'Na',      Na);
% setappdata(fig,'proc',    proc);
% setappdata(fig,'ax',      ax);
%
% drawTemperature(fig);
% set(fig,'KeyPressFcn',@temperatureCycleCallback);

% ╔══════════════════════════════════════════════════════════════════════╗
% ║ INTERACTIVE TEMPERATURE BROWSER                                     ║
% ╚══════════════════════════════════════════════════════════════════════╝

if config.steppedsine.system.BeamType == 2
    fig = figure('Name','Temperature Browser - Monolithic');
elseif config.steppedsine.system.BeamType == 3
    fig = figure('Name','Temperature Browser - Jointed');
else
    error('ERROR')
end
ax  = axes(fig);
hold(ax,'on'); grid(ax,'on');

xlabel(ax,'Frequency (Hz)');
ylabel(ax,'Response Amplitude');
set(ax,'YScale','log');
xlim([config.steppedsine.freq.Start config.steppedsine.freq.End]); ylim([0 100]);

nT    = numel(config.steppedsine.temperature.Vector);

% Storage for plot handles
hAmp  = gobjects(Na,1);
hBB   = gobjects(1);

activeZ = 1;

% Store shared data for callback
setappdata(fig,'activeZ',activeZ);
setappdata(fig,'temps',config.steppedsine.temperature.Vector);
setappdata(fig,'Na',Na);
setappdata(fig,'proc',proc);
setappdata(fig,'ax',ax);

% ============================================================
% SECOND FIGURE — FRF BROWSER (independent window)
% ============================================================

figFRF = figure('Name','FRF Browser');
axFRF  = axes(figFRF);
hold(axFRF,'on'); grid(axFRF,'on');

xlabel(axFRF,'Frequency (Hz)');
ylabel(axFRF,'|H(f)|');
set(axFRF,'YScale','log');

xlim(axFRF,[config.steppedsine.freq.Start config.steppedsine.freq.End]);

% Store required data
setappdata(figFRF,'proc',proc);
setappdata(figFRF,'temps',config.steppedsine.temperature.Vector);
setappdata(figFRF,'Na',Na);
setappdata(figFRF,'ax',axFRF);

% Link both figures together
setappdata(fig,'figFRF',figFRF);
setappdata(figFRF,'figMain',fig);

% Initial draw
drawTemperature(fig);

% --- Update FRF browser if it exists ---
if isappdata(fig,'figFRF')
    figFRF = getappdata(fig,'figFRF');
elseif isappdata(fig,'figMain')
    figFRF = fig;
    fig    = getappdata(fig,'figMain');
else
    return
end

if isvalid(figFRF)
    drawFRF(figFRF);
end

% Key control
set(fig,'KeyPressFcn',@temperatureCycleCallback);
set(figFRF,'KeyPressFcn',@temperatureCycleCallback);

% =========================================================================
% INDIVIDUAL BACKBONE PLOTS — one per temperature, using browser's draw logic
% =========================================================================

bbDir = fullfile('exports', 'figures', 'backbones', systemLabel);
if ~exist(bbDir, 'dir'); mkdir(bbDir); end

for zi = 1:nT

    T_val = config.steppedsine.temperature.Vector(zi);
    bb    = proc.steppedsine.backbone{zi};

    figBB = figure('Name', sprintf('%s  T = %d°C', systemLabel, T_val), ...
                   'Visible', 'off');
    ax_bb = axes(figBB);
    hold(ax_bb, 'on');
    grid(ax_bb, 'on');
    set(ax_bb, 'YScale', 'log');
    xlim(ax_bb, sort([config.steppedsine.freq.Start, config.steppedsine.freq.End]));
    ylim(ax_bb, [0 100]);
    xlabel(ax_bb, 'Frequency (Hz)');
    ylabel(ax_bb, 'Response Amplitude');
    % title(ax_bb, sprintf('%s — T = %d°C', systemLabel, T_val));

    % Replicate drawTemperature: plot each amplitude's peak-fit curve
    for a = 1:Na
        amp_a = proc.steppedsine.amp{zi}(a);
        if ~isempty(amp_a.Fit.fDense) && ~isempty(amp_a.Fit.yFit)
            plot(ax_bb, amp_a.Fit.fDense, amp_a.Fit.yFit, 'b-', 'LineWidth', 1);
        end
    end

    % % Backbone points
    % validBB = ~isnan(bb.f_pts) & ~isnan(bb.y_pts);
    % if any(validBB)
    %     scatter(ax_bb, bb.f_pts(validBB), bb.y_pts(validBB), 40, 'r', ...
    %             'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    % end

    % Backbone curve
    if ~all(isnan(bb.f)) && ~all(isnan(bb.y))
        plot(ax_bb, bb.f, bb.y, 'r-', 'LineWidth', 2);
    end

    exportgraphics(figBB, ...
        fullfile(bbDir, sprintf('backbone_%s_T%+04d.png', systemLabel, T_val)), ...
        'Resolution', 150);
    % close(figBB);

end
fprintf('Backbone figures saved to %s\n', bbDir);

fprintf('Backbone figures saved to %s\n', bbDir);

if config.steppedsine.system.BeamType == 2
    proc.backbone.monolithic = proc.steppedsine.backbone;
else
    proc.backbone.jointed = proc.steppedsine.backbone;
end
%% ==== END SteppedSine2.m ====

proc_all.Monolithic.steppedsine = proc.steppedsine;
proc_all.Monolithic.thermo      = proc.thermo;
% hysteresis.Monolithic.upsweep	= proc.steppedsine;
% config.steppedsine.data.fileName = arrayfun(@(k) ...
%     sprintf("0001-00%02d_152-147-0.1_1000_2000_2.csv", k), 1:config.steppedsine.nTemps);
% run SteppedSine2.m
% hysteresis.Monolithic.downsweep	= proc.steppedsine;

config.system.label    = 'Jointed';
config.system.BeamType = 3;
config.steppedsine.data.fileName = arrayfun(@(k) ...
    sprintf("0000-00%02d_112-117-0.1_1000_2000_3.csv", k), 1:config.steppedsine.nTemps);

%% ==== BEGIN SteppedSine2.m (integrated) ====

clear c1_fit  % Preserve exportCache across systems for combined results.
config.steppedsine.temperature.Vector = -30:10:120;
proc.steppedsine.backbone = cell(config.steppedsine.nTemps, 1);
proc.steppedsine.amp      = cell(config.steppedsine.nTemps, 1);
proc.steppedsine.FRF      = cell(config.steppedsine.nTemps, 1);

fJ0 = proc.impact.jointed.peakFreq;
fM0 = proc.impact.monolithic.peakFreq;

% Pre-loop preallocation
nTemps_pre = numel(-30:10:120);
Na_pre     = numel(0.5:0.5:8);

AMatrix        = NaN(Na_pre, nTemps_pre);
kEffMatrix     = NaN(Na_pre, nTemps_pre);
cEqMatrix      = NaN(Na_pre, nTemps_pre);
zetaMatrix     = NaN(Na_pre, nTemps_pre);
peakFreqMatrix = NaN(Na_pre, nTemps_pre);
peakAmpMatrix  = NaN(Na_pre, nTemps_pre);

summaryNonlinearityType = cell(nTemps_pre, 1);
summaryOnset            = NaN(1, nTemps_pre);
summaryOnsetUm          = NaN(1, nTemps_pre);

% Legacy descriptive onset threshold
isMonolithic_early = strcmp(config.system.label, 'Monolithic');
if isMonolithic_early
    onset_thresh_pct = 0.20;   % legacy descriptive threshold for monolithic case
else
    onset_thresh_pct = 0.50;   % legacy descriptive threshold for jointed case
end

% =========================================================================
% TEMPERATURE LOOP
% =========================================================================

for Z = 1:config.steppedsine.nTemps

    config.steppedsine.data.Name = config.steppedsine.data.fileName(Z);
    config.steppedsine.data.File = fullfile('data', config.steppedsine.data.Name);
    fprintf('Reading %s...\n', config.steppedsine.data.Name);
    config.steppedsine.data.opts = detectImportOptions(config.steppedsine.data.File);
    config.steppedsine.data.opts.DataLines = [24, Inf];
    config.steppedsine.data.Table = readtable(config.steppedsine.data.File, config.steppedsine.data.opts);

    % Parse metadata from filename
    p = sscanf(config.steppedsine.data.Name, '%d-%d_%d-%d-%f_%d_%d_%d');
    if numel(p) < 8; error('Filename parsing failed.'); end

    config.steppedsine.freq.Start               = p(3);
    config.steppedsine.freq.End                 = p(4);
    config.steppedsine.freq.Step                = p(5);
    config.steppedsine.sampling.Rate            = p(6);
    config.steppedsine.sampling.SamplesPerBlock = p(7);
    config.steppedsine.system.BeamType          = p(8);
    config.steppedsine.sampling.Transient       = 0;
    config.steppedsine.rows.N                   = height(config.steppedsine.data.Table);

    config.steppedsine.NumOfBlocks = floor(config.steppedsine.rows.N / config.steppedsine.sampling.SamplesPerBlock);

    % Amplitude settings
    config.steppedsine.amplitude.Min    = 0.5;
    config.steppedsine.amplitude.Max    = 8;
    config.steppedsine.amplitude.Step   = 0.5;
    config.steppedsine.amplitude.Vector = 0.5:0.5:8;
    config.steppedsine.amplitude.N      = numel(config.steppedsine.amplitude.Vector);
    Na = config.steppedsine.amplitude.N;
    config.steppedsine.voltage.Vector = linspace(8/Na, 8, Na);

    if mod(config.steppedsine.rows.N, Na) ~= 0
        error('Total rows not divisible by number of amplitude blocks.');
    end
    config.steppedsine.rows.PerAmp = config.steppedsine.rows.N / Na;
    config.steppedsine.freq.N      = round(config.steppedsine.rows.PerAmp / config.steppedsine.sampling.SamplesPerBlock);

    % System-dependent settings
    switch config.steppedsine.system.BeamType
        case 1; config.steppedsine.system.Mass = 0.71345;
        case 2; config.steppedsine.system.Mass = 0.22821;
        case 3; config.steppedsine.system.Mass = 0.23304;
        otherwise; error('Incorrect beam type.');
    end
    config.steppedsine.channels.ResponseStart = 3;
    config.steppedsine.channels.ResponseEnd   = 6;
    config.steppedsine.channels.Excitation    = 2;
    config.steppedsine.channels.Response      = 3:6;

    const.g  = 9.80665;
    const.fc = 35;

    Nf = config.steppedsine.freq.N;

    % Per-temperature preallocation
    proc.FRF.X    = NaN(Na, Nf);
    proc.FRF.Y    = NaN(Na, Nf);
    proc.FRF.H    = NaN(Na, Nf);
    proc.FRF.aMag = NaN(Na, Nf);
    proc.FRF.dMag = NaN(Na, Nf);
    proc.FRF.HMag    = NaN(Na, Nf);
    proc.FRF.HMagS   = NaN(Na, Nf);
    proc.FRF.HPhase  = NaN(Na, Nf);
    proc.FRF.HPhaseS = NaN(Na, Nf);

    proc.peaks.ResponseAcc = NaN(Na, Nf);
    proc.energy.Dissipated = NaN(Na, Nf);

    proc.amp = repmat(struct( ...
        'Voltage',       [], ...
        'PeakAmp',       NaN, ...
        'Zeta',          NaN, ...
        'DampingMethod', '', ...
        'HPB',           struct('Low',NaN,'High',NaN), ...
        'Circle',        struct('centre',NaN,'radius',NaN,'quality',NaN,'nPoints',0), ...
        'Lean',          struct('Left',NaN,'Right',NaN,'Asym',NaN), ...
        'Energy',        struct('Mean',NaN,'EDR',NaN,'NDI',NaN), ...
        'Fit',           struct('fDense',[],'yFit',[]), ...
        'Stiffness',     struct('k_eff',NaN,'k_ref',NaN,'deltaK',NaN)), Na, 1);

    backbone.response.f = NaN(1, Na);
    backbone.response.y = NaN(1, Na);

    proc.stiffness.f_res  = NaN(Na, 1);
    proc.stiffness.k_eff  = NaN(Na, 1);
    proc.stiffness.deltaK = NaN(Na, 1);

    fTest = (config.steppedsine.freq.Start : config.steppedsine.freq.Step : ...
             config.steppedsine.freq.Start + (Nf-1)*config.steppedsine.freq.Step)';

    % Detect descending sweep, store direction, then normalise Start/End
    sweepDir = sign(mean(diff( ...
        config.steppedsine.data.Table{1:config.steppedsine.sampling.SamplesPerBlock:end,1})));
    if sweepDir < 0
        blockOrder = reshape(1:config.steppedsine.rows.N, config.steppedsine.sampling.SamplesPerBlock, []);
        flippedIdx = blockOrder(:, end:-1:1); flippedIdx = flippedIdx(:);
        config.steppedsine.data.Table = config.steppedsine.data.Table(flippedIdx, :);
        fTest = flipud(fTest);
        % Normalise so Start < End (fixes xlim and any other range-dependent code)
        freqLims = sort([config.steppedsine.freq.Start, config.steppedsine.freq.End]);
        config.steppedsine.freq.Start = freqLims(1);
        config.steppedsine.freq.End   = freqLims(2);
    end
    % Store sweep direction for hysteresis analysis (+1 = up, -1 = down)
    config.steppedsine.sweep.direction = sweepDir;

    % Retrieve impact-test reference values.
    % k_ref_impact comes from the system-specific impact test (correct — each
    % system's own baseline stiffness at T=20°C).
    % m_eff uses the MONOLITHIC physical mass for both systems so that k_eff,
    % c_eq and all regression coefficients are on the same absolute scale and
    % the jointed system correctly appears softer than the monolithic.
    if strcmp(config.system.label, 'Monolithic')
        k_ref_impact = proc.impact.monolithic.stiffness;
    else
        k_ref_impact = proc.impact.jointed.stiffness;
    end
    m_eff = proc.impact.monolithic.m_eff;   % shared physical mass — same for both systems
    T_impact  = 20;
    T_current = config.steppedsine.temperature.Vector(Z);

% =========================================================================
% AMPLITUDE LOOP
% =========================================================================

    for a = 1:Na

        proc.amp(a).Voltage = config.steppedsine.voltage.Vector(a);
        aStart = (a-1) * config.steppedsine.rows.PerAmp + 1;

        % Frequency loop — signal processing and FRF
        for k = 1:Nf
            r0 = aStart + (k-1)*config.steppedsine.sampling.SamplesPerBlock;
            r1 = r0 + config.steppedsine.sampling.SamplesPerBlock - 1;

            t   = config.steppedsine.data.Table{r0:r1, 1};
            t   = t - t(1);
            exc = detrend(config.steppedsine.data.Table{r0:r1, config.steppedsine.channels.Excitation});
            res = detrend(config.steppedsine.data.Table{r0:r1, config.steppedsine.channels.Response});
            res = mean(res, 2);

            fs = 1 / mean(diff(t));
            [bB, aB] = butter(2, const.fc/(fs/2), 'high');
            res = filtfilt(bB, aB, res);

            nt  = round(config.steppedsine.sampling.Transient * config.steppedsine.sampling.SamplesPerBlock);
            t   = t(nt+1:end);
            exc = exc(nt+1:end);
            res = res(nt+1:end);

            vel  = cumtrapz(t, res);
            vel  = vel - polyval(polyfit(t,vel,1), t);
            dis  = cumtrapz(t, vel);
            dis  = dis - polyval(polyfit(t,dis,2), t);

            proc.energy.Dissipated(a,k) = trapz(t, exc .* vel);

            fExc  = fTest(k);
            omega = 2*pi*fExc;

            proc.FRF.X(a,k) = 2/numel(exc) * sum(exc .* exp(-1i*2*pi*fExc*t));
            proc.FRF.Y(a,k) = 2/numel(res) * sum(res .* exp(-1i*2*pi*fExc*t));
            % Pipeline-normalised complex response quantity. Because the denominator
            % includes the specimen-mass and g scaling used in the acquisition workflow,
            % this is not automatically interpreted as a calibrated force/acceleration FRF.
            proc.FRF.H(a,k) = proc.FRF.Y(a,k) / (config.steppedsine.system.Mass * const.g * proc.FRF.X(a,k));

            proc.FRF.aMag(a,k) = abs(proc.FRF.Y(a,k));
            proc.FRF.dMag(a,k) = proc.FRF.aMag(a,k) / omega^2;

            proc.peaks.ResponseAcc(a,k) = localPeakEstimate(res);
        end

        proc.FRF.HMag(a,:)    = abs(proc.FRF.H(a,:));
        proc.FRF.HMagS(a,:)   = smoothdata(proc.FRF.HMag(a,:),  'movmean', 3);
        proc.FRF.HPhase(a,:)  = angle(proc.FRF.H(a,:));
        proc.FRF.HPhaseS(a,:) = smoothdata(proc.FRF.HPhase(a,:),'movmean', 3);

        % Resonance peak — PCHIP on acceleration response
        y = proc.peaks.ResponseAcc(a,:);
        proc.amp(a).Fit.fRaw = fTest(:);
        proc.amp(a).Fit.yRaw = y(:);
        if all(isnan(y)) || all(y == 0); continue; end

        validPts = ~isnan(y);
        pp = pchip(fTest(validPts), y(validPts));
        fD = linspace(min(fTest(validPts)), max(fTest(validPts)), ...
                      max(200, round((max(fTest(validPts))-min(fTest(validPts)))/0.001)));
        yF = ppval(pp, fD);

        [yPk, iPk] = max(yF);
        fPk        = fD(iPk);
        omega_pk   = 2*pi*fPk;

        proc.amp(a).PeakAmp   = yPk;
        peakAmpMatrix(a, Z)   = yPk;
        peakFreqMatrix(a, Z)  = fPk;
        proc.amp(a).Fit.fDense = fD;
        proc.amp(a).Fit.yFit   = yF;
        proc.amp(a).Energy.Mean = mean(abs(proc.energy.Dissipated(a,:)));

        % Effective stiffness
        if exist('c1_fit','var') && ~isnan(c1_fit)
            k_ref_T = k_ref_impact * (1 + c1_fit * (T_current - T_impact));
        else
            k_ref_T = k_ref_impact;
        end
        k_eff_loc = m_eff * omega_pk^2;

        proc.amp(a).Stiffness.k_eff  = k_eff_loc;
        proc.amp(a).Stiffness.k_ref  = k_ref_T;
        proc.amp(a).Stiffness.deltaK = (k_eff_loc - k_ref_T) / k_ref_T;
        proc.stiffness.f_res(a)  = fPk;
        proc.stiffness.k_eff(a)  = k_eff_loc;
        proc.stiffness.deltaK(a) = proc.amp(a).Stiffness.deltaK;

        % Harmonic displacement estimate reconstructed from acceleration at fPk
        HdMag = proc.FRF.dMag(a, :);
        validH_d = ~isnan(HdMag);
        if sum(validH_d) >= 3
            dMag_atPk = ppval(pchip(fTest(validH_d), HdMag(validH_d)), fPk);
        else
            [~, iPk_snap] = min(abs(fTest - fPk));
            dMag_atPk = HdMag(iPk_snap);
        end
        A_metres = dMag_atPk * const.g;
        fprintf('A_metres (a=%d, Z=%d): %.6e m  (%.4f um)\n', a, Z, A_metres, A_metres*1e6);
        AMatrix(a, Z)    = A_metres;
        kEffMatrix(a, Z) = k_eff_loc;

        % Damping — Nyquist circle fit, fallback to half-power bandwidth
        Hc = proc.FRF.H(a, :);
        HmagFull = abs(Hc);
        validHm  = ~isnan(HmagFull) & HmagFull > 0;
        if sum(validHm) >= 3
            HmFit  = ppval(pchip(fTest(validHm), HmagFull(validHm)), fD);
            aboveH = HmFit >= max(HmFit)/sqrt(2);
            if any(aboveH)
                bw_est = fD(find(aboveH,1,'last')) - fD(find(aboveH,1,'first'));
                bw_est = max(bw_est, 4*config.steppedsine.freq.Step);
            else
                bw_est = 4*config.steppedsine.freq.Step;
            end
        else
            bw_est = 4*config.steppedsine.freq.Step;
        end

        bw_factor = 2.5;
        win_mask  = fTest' >= (fPk - bw_factor*bw_est) & ...
                    fTest' <= (fPk + bw_factor*bw_est) & ...
                    ~isnan(Hc) & abs(Hc) > 0;
        Hc_win = Hc(win_mask);
        f_win  = fTest(win_mask);

        zeta_circle = NaN;
        circle_ok   = false;

        if numel(Hc_win) >= 5
            xp = real(Hc_win(:));
            yp = imag(Hc_win(:));
            coef = [xp, yp, ones(size(xp))] \ -(xp.^2 + yp.^2);
            xc = -coef(1)/2;   yc = -coef(2)/2;
            R2c = xc^2 + yc^2 - coef(3);
            if R2c > 0
                Rc = sqrt(R2c);
                dist_from_circle = abs(sqrt((xp-xc).^2 + (yp-yc).^2) - Rc);
                circle_quality   = mean(dist_from_circle) / (Rc + eps);

                [~, i_near] = min(abs(f_win - fPk));
                i_lo2 = max(1, i_near-1);
                i_hi2 = min(numel(f_win), i_near+1);

                if i_lo2 ~= i_hi2 && circle_quality < 0.3
                    dtheta = abs(angle(Hc_win(i_hi2) - complex(xc,yc)) - ...
                                 angle(Hc_win(i_lo2) - complex(xc,yc)));
                    if dtheta > pi; dtheta = 2*pi - dtheta; end
                    domega = 2*pi*(f_win(i_hi2) - f_win(i_lo2));
                    if domega > 0 && dtheta > 0
                        zeta_circle = 1 / (omega_pk * (dtheta/domega));
                        circle_ok = zeta_circle > 0.001 && zeta_circle < 0.15;
                    end
                end

                proc.amp(a).Circle.centre  = complex(xc, yc);
                proc.amp(a).Circle.radius  = Rc;
                proc.amp(a).Circle.quality = circle_quality;
                proc.amp(a).Circle.nPoints = numel(Hc_win);
            end
        end

        if circle_ok
            proc.amp(a).Zeta          = zeta_circle;
            proc.amp(a).DampingMethod = 'CircleFit';
        else
            yH = yPk / sqrt(2);
            f1 = findHalfPowerIntersection(fD, yF, iPk, 'left',  yH);
            f2 = findHalfPowerIntersection(fD, yF, iPk, 'right', yH);
            if ~isnan(f1) && ~isnan(f2) && f2 > f1
                proc.amp(a).Zeta          = (f2 - f1) / (2*fPk);
                proc.amp(a).HPB.Low       = f1;
                proc.amp(a).HPB.High      = f2;
                proc.amp(a).DampingMethod = 'HalfPower';
            else
                proc.amp(a).Zeta          = NaN;
                proc.amp(a).DampingMethod = 'Failed';
            end
        end

        % Peak asymmetry (always from PCHIP)
        yH = yPk / sqrt(2);
        f1l = findHalfPowerIntersection(fD, yF, iPk, 'left',  yH);
        f2l = findHalfPowerIntersection(fD, yF, iPk, 'right', yH);
        if ~isnan(f1l) && ~isnan(f2l) && f2l > f1l
            proc.amp(a).Lean.Left  = (fPk - f1l) / (f2l - f1l);
            proc.amp(a).Lean.Right = (f2l - fPk) / (f2l - f1l);
            proc.amp(a).Lean.Asym  = (f2l + f1l - 2*fPk) / (f2l - f1l);
        end

        backbone.response.f(a) = fPk;
        backbone.response.y(a) = yPk;

        zetaMatrix(a, Z) = proc.amp(a).Zeta;
        if ~isnan(proc.amp(a).Zeta)
            cEqMatrix(a, Z) = 2 * m_eff * proc.amp(a).Zeta * omega_pk;
        end

    end  % amplitude loop

    % =========================================================================
% EXTRACT EXCITATION ACCELERATION LEVELS (g RMS and peak)
% =========================================================================

exc_g_rms  = NaN(Na, 1);
exc_g_peak = NaN(Na, 1);

for a = 1:Na
    aStart = (a-1) * config.steppedsine.rows.PerAmp + 1;
    aEnd   = aStart + config.steppedsine.rows.PerAmp - 1;

    % Pull the full excitation block for this amplitude (all frequencies)
    exc_raw = config.steppedsine.data.Table{aStart:aEnd, ...
                  config.steppedsine.channels.Excitation};
    exc_raw = detrend(exc_raw);

    % Convert from raw units to g
    % (your data is in m/s^2 based on const.g usage — divide to get g)
    exc_g = exc_raw;% / const.g;

    exc_g_rms(a)  = rms(exc_g);
    exc_g_peak(a) = max(abs(exc_g));
end

fprintf('\nExcitation levels at T = %d°C:\n', T_current);
fprintf('  Amp level | RMS (g) | Peak (g)\n');
for a = 1:Na
    fprintf('  %9.1f | %7.4f | %8.4f\n', ...
        config.steppedsine.amplitude.Vector(a), exc_g_rms(a), exc_g_peak(a));
end

    fprintf('...done.\n\n');

    % Legacy uncalibrated energy-proxy normalisation (not a physical dissipation measurement)
    E0 = proc.amp(1).Energy.Mean;
    for a = 1:Na
        proc.amp(a).Energy.EDR = proc.amp(a).Energy.Mean / (E0 + eps);
        if ~isnan(proc.amp(a).PeakAmp) && proc.amp(a).PeakAmp > 0
            proc.amp(a).Energy.NDI = proc.amp(a).Energy.Mean / proc.amp(a).PeakAmp^2;
        end
        proc.amp(a).Voltage = config.steppedsine.voltage.Vector(a);
    end

    % Backbone polynomial fit (quadratic max)
    validBB = ~isnan(backbone.response.f) & ~isnan(backbone.response.y);
    if sum(validBB) >= 2
        aIdx   = find(validBB);
        pOrder = min(2, numel(aIdx)-1);
        p_fx   = polyfit(aIdx, backbone.response.f(aIdx), pOrder);
        p_fy   = polyfit(aIdx, backbone.response.y(aIdx), pOrder);
        aDense = linspace(1, max(aIdx), 10000);
        fCurve = polyval(p_fx, aDense);
        yCurve = polyval(p_fy, aDense);
        if sum(diff(sign(diff(fCurve))) ~= 0) > 1
            fCurve = polyval(polyfit(aIdx, backbone.response.f(aIdx), 1), aDense);
        end
    else
        fCurve = NaN;  yCurve = NaN;
    end

    proc.steppedsine.backbone{Z} = struct( ...
        'f',     fCurve(:), ...
        'y',     yCurve(:), ...
        'f_pts', backbone.response.f(:), ...
        'y_pts', backbone.response.y(:), ...
        'A_pts', AMatrix(:, Z), ...
        'k_pts', kEffMatrix(:, Z));
    proc.steppedsine.amp{Z} = proc.amp;
    proc.steppedsine.FRF{Z} = struct('HMagS', proc.FRF.HMagS, ...
                                      'HMag',  proc.FRF.HMag, ...
                                      'f',     fTest);

    % Legacy descriptive onset estimate
    dK     = proc.stiffness.deltaK;
    valid  = ~isnan(dK);
    onsetIdx = NaN;
    if sum(valid) >= 2
        exceed = abs(dK)*100 > onset_thresh_pct & valid;
        for a = 1:Na-1
            if exceed(a) && exceed(a+1)
                onsetIdx = a; break;
            end
        end
    end
    if isnan(onsetIdx)
        summaryOnset(Z)   = config.steppedsine.amplitude.Vector(Na);
        summaryOnsetUm(Z) = AMatrix(Na, Z) * 1e6;
    else
        summaryOnset(Z)   = config.steppedsine.amplitude.Vector(onsetIdx);
        summaryOnsetUm(Z) = AMatrix(onsetIdx, Z) * 1e6;
    end

    % Provisional descriptive label; replaced later by BFS classification
    dK_valid = dK(valid);
    summaryNonlinearityType{Z} = 'Softening';
    if mean(dK_valid) > 0; summaryNonlinearityType{Z} = 'Hardening'; end

end  % temperature loop

% =========================================================================
% DUFFING COEFFICIENT IDENTIFICATION
% f_res^2 = f0^2 + alpha * A^2  (cubic Duffing backbone linearisation)
% alpha > 0 = hardening;  alpha < 0 = softening
% =========================================================================

T_vec = config.steppedsine.temperature.Vector(:);
nT = config.steppedsine.nTemps;
duffing_alpha = NaN(1, nT);
duffing_f0sq  = NaN(1, nT);
duffing_R2    = NaN(1, nT);

for Z = 1:nT
    f_pts = proc.steppedsine.backbone{Z}.f_pts(:);   % resonance freq per amplitude
    A_pts = proc.steppedsine.backbone{Z}.A_pts(:);   % displacement per amplitude [m]
    valid = ~isnan(f_pts) & ~isnan(A_pts) & A_pts > 0;
    if sum(valid) < 4; continue; end

    f2 = f_pts(valid).^2;
    A2 = (A_pts(valid)*1e6).^2;   % convert to µm^2 for numerically stable fit

    X    = [ones(sum(valid),1), A2];
    beta = X \ f2;
    duffing_f0sq(Z)  = beta(1);
    duffing_alpha(Z) = beta(2);   % Hz^2 / µm^2
    fprintf('Duffing fit Z=%d: alpha=%.6f Hz2/um2  f0=%.4f Hz\n', Z, beta(2), sqrt(beta(1)));

    f2_pred = X * beta;
    ss_res  = sum((f2 - f2_pred).^2);
    ss_tot  = sum((f2 - mean(f2)).^2);
    if ss_tot > 0; duffing_R2(Z) = 1 - ss_res/ss_tot; end
end

proc.thermo.duffing.alpha  = duffing_alpha;
proc.thermo.duffing.f0sq   = duffing_f0sq;
proc.thermo.duffing.R2     = duffing_R2;
proc.thermo.duffing.T_vec  = T_vec(:)';

fprintf('Duffing alpha(T): min=%.3f  max=%.3f  Hz^2/um^2\n', ...
    min(duffing_alpha,[],'omitnan'), max(duffing_alpha,[],'omitnan'));
fprintf('R2 range: %.3f to %.3f\n', ...
    min(duffing_R2,[],'omitnan'), max(duffing_R2,[],'omitnan'));

% =========================================================================
% ALPHA(T) POLYNOMIAL FIT — Global Thermal-Mechanical Duffing Model
% f_res^2(A,T) = f0^2(T) + alpha(T) * A^2
% where alpha(T) = sum_k p_k * T^k  (polynomial in absolute temperature)
% and   f0^2(T)  = sum_k q_k * T^k
% =========================================================================

T_abs   = proc.thermo.duffing.T_vec(:);       % absolute temperature [°C]
alpha_v = proc.thermo.duffing.alpha(:);        % Hz^2/µm^2 per temperature
f0sq_v  = proc.thermo.duffing.f0sq(:);        % Hz^2 per temperature
R2_v    = proc.thermo.duffing.R2(:);

% --- Degree selection: fit poly orders 1–4, pick by adjusted R² ---
% Weight by per-temperature Duffing R² so unreliable alpha estimates
% (cold jointed, R²~0.03) contribute less to the polynomial fit.
w = max(R2_v, 0);   % zero-floor negative R² values

alpha_poly = struct();
f0sq_poly  = struct();

for deg = 1:4
    % Alpha(T) polynomial
    % Only positive-weight, finite observations contribute to weighted fitting.
    valid_a = isfinite(alpha_v) & isfinite(T_abs) & isfinite(w) & w > 0;
    n_alpha = nnz(valid_a);
    n_parameters = deg + 1;  % degree d polynomial includes the intercept
    if n_alpha <= n_parameters
        error('research_pipeline:InsufficientDuffingData', ...
            'Degree %d alpha(T) model needs >%d positive-weight observations; found %d.', ...
            deg, n_parameters, n_alpha);
    end
    p_a = wpolyfit(T_abs(valid_a), alpha_v(valid_a), w(valid_a), deg);
    alpha_pred = polyval(p_a, T_abs(valid_a));
    w_alpha = w(valid_a);
    y_alpha = alpha_v(valid_a);
    mean_alpha = sum(w_alpha .* y_alpha) / sum(w_alpha);
    ss_res = sum(w_alpha .* (y_alpha-alpha_pred).^2);
    ss_tot = sum(w_alpha .* (y_alpha-mean_alpha).^2);
    R2_a = 1 - ss_res / max(ss_tot, eps);
    % Adjusted R2 denominator is n - p, where p includes the intercept.
    adjR2_a = 1 - (1-R2_a)*(n_alpha-1)/(n_alpha-n_parameters);

    % f0sq(T) uses its own valid sample count and unweighted fit.
    valid_f = isfinite(f0sq_v) & isfinite(T_abs);
    n_f0 = nnz(valid_f);
    if n_f0 <= n_parameters
        error('research_pipeline:InsufficientDuffingData', ...
            'Degree %d f0sq(T) model needs >%d observations; found %d.', ...
            deg, n_parameters, n_f0);
    end
    p_f = polyfit(T_abs(valid_f), f0sq_v(valid_f), deg);
    f0sq_pred = polyval(p_f, T_abs(valid_f));
    ss_res_f = sum((f0sq_v(valid_f)-f0sq_pred).^2);
    ss_tot_f = sum((f0sq_v(valid_f)-mean(f0sq_v(valid_f))).^2);
    R2_f = 1 - ss_res_f / max(ss_tot_f, eps);
    adjR2_f = 1 - (1-R2_f)*(n_f0-1)/(n_f0-n_parameters);

    alpha_poly(deg).coeffs = p_a;
    alpha_poly(deg).R2     = R2_a;
    alpha_poly(deg).adjR2  = adjR2_a;
    f0sq_poly(deg).coeffs  = p_f;
    f0sq_poly(deg).R2      = R2_f;
    f0sq_poly(deg).adjR2   = adjR2_f;

    fprintf('deg %d | alpha(T): R²=%.4f adjR²=%.4f | f0²(T): R²=%.4f adjR²=%.4f\n', ...
        deg, R2_a, adjR2_a, R2_f, adjR2_f);
end

% --- Select best degree by adjusted R² for alpha(T) ---
adjR2_vals = [alpha_poly.adjR2];
[~, best_deg] = max(adjR2_vals);
p_alpha_best  = alpha_poly(best_deg).coeffs;
p_f0sq_best   = f0sq_poly(best_deg).coeffs;
fprintf('\nSelected polynomial degree: %d  (alpha adjR²=%.4f)\n', ...
    best_deg, alpha_poly(best_deg).adjR2);

% --- Global surface prediction: f_res^2(A,T) ---
% Predict over a dense (T, A) grid for plotting and export
T_dense = linspace(min(T_abs), max(T_abs), 200);
A_dense = linspace(0, max(AMatrix(:), [], 'omitnan')*1e6, 100);   % µm
[T_grid, A_grid] = meshgrid(T_dense, A_dense);

alpha_grid = polyval(p_alpha_best, T_grid);
f0sq_grid  = polyval(p_f0sq_best,  T_grid);
fres_grid  = sqrt(max(f0sq_grid + alpha_grid .* A_grid.^2, 0));

% --- Store results ---
proc.thermo.duffing.alpha_poly = alpha_poly;
proc.thermo.duffing.f0sq_poly  = f0sq_poly;
proc.thermo.duffing.best_deg   = best_deg;
proc.thermo.duffing.p_alpha    = p_alpha_best;
proc.thermo.duffing.p_f0sq     = p_f0sq_best;
proc.thermo.duffing.T_grid     = T_grid;
proc.thermo.duffing.A_grid     = A_grid;
proc.thermo.duffing.fres_grid  = fres_grid;

% --- Diagnostic figure ---
figure('Name', sprintf('Alpha(T) Polynomial — %s', config.system.label));

subplot(1,3,1);
scatter(T_abs(valid_a), alpha_v(valid_a), 60, w(valid_a), 'filled'); hold on;
T_plot = linspace(min(T_abs), max(T_abs), 300);
plot(T_plot, polyval(p_alpha_best, T_plot), 'r-', 'LineWidth', 2);
yline(0, 'k--');
xlabel('Temperature (°C)'); ylabel('\alpha (Hz^2/\mum^2)');
title(sprintf('\\alpha(T) — degree %d  adjR^2=%.3f', best_deg, alpha_poly(best_deg).adjR2));
colorbar; clabel_text = 'Duffing R^2'; % colour = per-T reliability
grid on;

subplot(1,3,2);
scatter(T_abs(valid_f), f0sq_v(valid_f), 60, 'k', 'filled'); hold on;
plot(T_plot, polyval(p_f0sq_best, T_plot), 'r-', 'LineWidth', 2);
xlabel('Temperature (°C)'); ylabel('f_0^2 (Hz^2)');
title(sprintf('f_0^2(T) — degree %d  R^2=%.3f', best_deg, f0sq_poly(best_deg).R2));
grid on;

subplot(1,3,3);
surf(T_grid, A_grid, fres_grid, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
xlabel('Temperature (°C)'); ylabel('Amplitude (\mum)'); zlabel('f_{res} (Hz)');
title('Global Model: f_{res}(A,T)');
colorbar; view(45, 30); grid on;

% =========================================================================
% MINDLIN POWER-LAW FIT: delta_K ~ C * A^n  above the zero-crossing
% Uses Duffing f0 as reference, consistent with Table 2 deltaK definition
% =========================================================================

powerlaw_n   = NaN(1, nT);
powerlaw_C   = NaN(1, nT);
powerlaw_R2  = NaN(1, nT);

for Z = 1:nT
    f_pts = proc.steppedsine.backbone{Z}.f_pts(:);
    A_pts = proc.steppedsine.backbone{Z}.A_pts(:) * 1e6;  % µm

    valid = ~isnan(f_pts) & ~isnan(A_pts) & A_pts > 0;
    if sum(valid) < 4; continue; end

    % Use Duffing-extrapolated f0 as reference — same as Table 2 deltaK
    f0_sq = duffing_f0sq(Z);
    if isnan(f0_sq) || f0_sq <= 0; continue; end
    f0_T  = sqrt(f0_sq);

    dK_vec = ((f_pts(valid).^2 - f0_sq) / f0_sq) * 100;  % % deviation
    A_vec  = A_pts(valid);
    fprintf('dK Z=%d: A range %.4f-%.4f um  dK range %.3f-%.3f%%\n', Z, min(A_vec), max(A_vec), min(dK_vec), max(dK_vec));

    % Find crossover: dK transitions from negative to positive
    cross_idx = NaN;
    for ai = 1:numel(dK_vec)-1
        if dK_vec(ai) < 0 && dK_vec(ai+1) >= 0
            cross_idx = ai + 1; break;
        end
    end
    if isnan(cross_idx); continue; end

    A_above  = A_vec(cross_idx:end);
    dK_above = dK_vec(cross_idx:end);
    pos      = dK_above > 0;
    if sum(pos) < 3; continue; end

    logA  = log(A_above(pos));
    logdK = log(dK_above(pos));
    p_fit = polyfit(logA, logdK, 1);
    n_fit = p_fit(1);
    C_fit = exp(p_fit(2));

    dK_pred = C_fit * A_above(pos).^n_fit;
    ss_res  = sum((dK_above(pos) - dK_pred).^2);
    ss_tot  = sum((dK_above(pos) - mean(dK_above(pos))).^2);

    powerlaw_n(Z)  = n_fit;
    powerlaw_C(Z)  = C_fit;
    powerlaw_R2(Z) = 1 - ss_res / (ss_tot + eps);
end

proc.thermo.powerlaw.n   = powerlaw_n;
proc.thermo.powerlaw.C   = powerlaw_C;
proc.thermo.powerlaw.R2  = powerlaw_R2;

% Print table
fprintf('\nMindlin power-law exponent n(T):\n');
fprintf('%-10s %-8s %-8s %-8s\n', 'T (degC)', 'n', 'C', 'R2');
for Z = 1:nT
    T_val = config.steppedsine.temperature.Vector(Z);
    if ~isnan(powerlaw_n(Z))
        fprintf('%-10d %-8.3f %-8.4f %-8.4f\n', ...
            T_val, powerlaw_n(Z), powerlaw_C(Z), powerlaw_R2(Z));
    else
        fprintf('%-10d  — (no crossover)\n', T_val);
    end
end
fprintf('Mean n = %.3f,  std = %.3f\n', ...
    mean(powerlaw_n,'omitnan'), std(powerlaw_n,'omitnan'));

% =========================================================================
% REGRESSION
% =========================================================================

T_vec    = config.steppedsine.temperature.Vector(:);
T0       = 20;
dT_vec   = T_vec - T0;
nT       = numel(T_vec);

A_flat    = AMatrix(:);
dT_flat   = repmat(dT_vec', Na, 1);  dT_flat  = dT_flat(:);
T_flat    = repmat(T_vec',  Na, 1);  T_flat   = T_flat(:);
kEff_flat = kEffMatrix(:);
cEq_flat  = cEqMatrix(:);

valid_k     = ~isnan(A_flat) & ~isnan(kEff_flat);
valid_c_raw = ~isnan(A_flat) & ~isnan(cEq_flat);

A_v  = A_flat(valid_k);   dT_v = dT_flat(valid_k);   k_v = kEff_flat(valid_k);
A_c_all  = A_flat(valid_c_raw);
dT_c_all = dT_flat(valid_c_raw);
T_c_all  = T_flat(valid_c_raw);
c_v_all  = cEq_flat(valid_c_raw);

% -------------------------------------------------------------------------
% STIFFNESS-RELATED REGRESSION — linear-T piecewise
% Model: k = k_ref * (1 + c1*dT + [g0_regime + g1_regime*dT]*A^2)
% c2*dT^2 removed: degraded R² in both systems.
%
% T* selection: use fixed value when config.steppedsine.T_star_fixed is provided
% is set; otherwise sweep 0:10:60 and pick best R².
% -------------------------------------------------------------------------

use_fixed_Tstar = isfield(config.steppedsine,'T_star_fixed') && ...
                  ~isnan(config.steppedsine.T_star_fixed);

T_star_candidates = 0:10:60;
R2_k_sweep        = NaN(size(T_star_candidates));

if use_fixed_Tstar
    T_star_opt = config.steppedsine.T_star_fixed - T0;
    fprintf('\n=== T* fixed at %.0f°C (dT=%.0f) — R² sweep skipped ===\n', ...
        config.steppedsine.T_star_fixed, T_star_opt);
    r_fixed = double(dT_v >= T_star_opt);
    if sum(r_fixed<1)>=10 && sum(r_fixed>0)>=10
        Phi_f = [ones(size(A_v)), dT_v, ...
                 (1-r_fixed).*A_v.^2, (1-r_fixed).*dT_v.*A_v.^2, ...
                 r_fixed.*A_v.^2,      r_fixed.*dT_v.*A_v.^2];
        k_hat_f = Phi_f * (Phi_f \ k_v);
        R2_k = 1 - sum((k_v-k_hat_f).^2)/sum((k_v-mean(k_v)).^2);
    else
        R2_k = NaN;
        warning('Insufficient data in one regime for fixed T*=%.0f°C.', ...
            config.steppedsine.T_star_fixed);
    end
else
    for iStar = 1:numel(T_star_candidates)
        T_s = T_star_candidates(iStar);
        r_s = double(dT_v >= T_s);
        if sum(r_s<1)<10 || sum(r_s>0)<10; continue; end
        Phi_s = [ones(size(A_v)), dT_v, ...
                 (1-r_s).*A_v.^2, (1-r_s).*dT_v.*A_v.^2, ...
                 r_s.*A_v.^2,      r_s.*dT_v.*A_v.^2];
        k_hat = Phi_s * (Phi_s \ k_v);
        R2_k_sweep(iStar) = 1 - sum((k_v-k_hat).^2)/sum((k_v-mean(k_v)).^2);
    end
    [R2_k, iStar_best] = max(R2_k_sweep);
    T_star_opt = T_star_candidates(iStar_best);
end

regime_k = double(dT_v >= T_star_opt);

% Final linear-T piecewise fit
Phi_k  = [ones(size(A_v)), dT_v, ...
          (1-regime_k).*A_v.^2, (1-regime_k).*dT_v.*A_v.^2, ...
          regime_k.*A_v.^2,      regime_k.*dT_v.*A_v.^2];
theta_k = Phi_k \ k_v;
k_pred  = Phi_k * theta_k;

k_ref_fit   = theta_k(1);
c1_fit      = theta_k(2) / k_ref_fit;
gamma0_cold = theta_k(3) / k_ref_fit;
gamma1_cold = theta_k(4) / k_ref_fit;
gamma0_warm = theta_k(5) / k_ref_fit;
gamma1_warm = theta_k(6) / k_ref_fit;

% Global model (no piecewise split)
Phi_k_g    = [ones(size(A_v)), dT_v, A_v.^2, dT_v.*A_v.^2];
theta_k_g  = Phi_k_g \ k_v;
R2_k_global = 1 - sum((k_v - Phi_k_g*theta_k_g).^2) / sum((k_v-mean(k_v)).^2);
gamma0_fit  = theta_k_g(3)/theta_k_g(1);
gamma1_fit  = theta_k_g(4)/theta_k_g(1);

% Stiffness-related residual structure by temperature and amplitude index
k_resid = k_v - k_pred;
T_resid_mean = NaN(1,nT);  T_resid_std = NaN(1,nT);
A_resid_mean = NaN(1,Na);  A_resid_std = NaN(1,Na);
for iT_r = 1:nT
    idx_r = dT_v == dT_vec(iT_r);
    if sum(idx_r)>=2
        T_resid_mean(iT_r) = mean(k_resid(idx_r));
        T_resid_std(iT_r)  = std(k_resid(idx_r));
    end
end
A_vec_flat_idx = repmat((1:Na)',nT,1);
A_v_idx        = A_vec_flat_idx(valid_k);
for iA_r = 1:Na
    idx_r = A_v_idx == iA_r;
    if sum(idx_r)>=2
        A_resid_mean(iA_r) = mean(k_resid(idx_r));
        A_resid_std(iA_r)  = std(k_resid(idx_r));
    end
end

if strcmp(config.system.label,'Monolithic')
    f_ref_Hz = fM0;
else
    f_ref_Hz = fJ0;
end

fprintf('\n=== STIFFNESS (piecewise, linear-T) ===\n');
if use_fixed_Tstar
    fprintf('T_star=%.0f°C (FIXED)  k_ref=%.2f N/m  f_ref=%.4f Hz\n', T0+T_star_opt, k_ref_fit, f_ref_Hz);
else
    fprintf('T_star=%.0f°C (swept)  k_ref=%.2f N/m  f_ref=%.4f Hz\n', T0+T_star_opt, k_ref_fit, f_ref_Hz);
end
fprintf('  c1=%.4e /°C\n', c1_fit);
f_ratio_sq = (fJ0/fM0)^2;
fprintf('Sanity: (f_J/f_M)^2=%.4f  joint=%.1f%% of mono\n', f_ratio_sq, 100*f_ratio_sq);
fprintf('Cold: g0=%.4e  g1=%.4e  /m²\n', gamma0_cold, gamma1_cold);
fprintf('Warm: g0=%.4e  g1=%.4e  /m²\n', gamma0_warm, gamma1_warm);
fprintf('R²: piecewise=%.4f  global=%.4f\n', R2_k, R2_k_global);
[~,iT_maxres] = max(abs(T_resid_mean),[],'omitnan');
fprintf('Max |residual mean| by T: %.1f N/m at T=%.0f°C\n', ...
    abs(T_resid_mean(iT_maxres)), T_vec(iT_maxres));

% -------------------------------------------------------------------------
% EQUIVALENT-DAMPING REGRESSION — bounded exponential model
% Model: c_eq = c_ref*exp(d1*dT) + d2*A + d3*A^2
% Two-stage: (1) log-scale regression for d1 from low-amplitude data;
% (2) amplitude residuals for d2, d3.
% Hard bound: |d1| <= d1_max = 0.05 /°C.
% -------------------------------------------------------------------------

MAD_thresh     = 3.0;
c_outlier_flag = false(size(c_v_all));
for iT = 1:numel(unique(T_c_all))
    T_lvl = unique(T_c_all); T_lvl = T_lvl(iT);
    idx_T = T_c_all == T_lvl;
    c_T   = c_v_all(idx_T);
    c_med = median(c_T);
    mad_T = median(abs(c_T - c_med));
    if mad_T > eps
        c_outlier_flag(idx_T) = abs(c_T - c_med) > MAD_thresh * mad_T;
    end
end
n_flagged = sum(c_outlier_flag);
fprintf('\n=== DAMPING (MAD outliers: %d removed of %d) ===\n', n_flagged, numel(c_v_all));

valid_c = ~c_outlier_flag;
A_c  = A_c_all(valid_c);
dT_c = dT_c_all(valid_c);
T_c  = T_c_all(valid_c);
c_v  = c_v_all(valid_c);

isMonolithic = strcmp(config.system.label, 'Monolithic');

% Baseline comparison (unchanged model, for R² tracking)
Phi_c_old  = [ones(size(A_c)), dT_c, A_c, dT_c.*A_c];
c_pred_old = Phi_c_old * (Phi_c_old \ c_v);
R2_c_old   = 1 - sum((c_v-c_pred_old).^2)/sum((c_v-mean(c_v)).^2);

regime_c = double(T_c >= T_star_opt + T0);
cold_c   = 1 - regime_c;
warm_c   =     regime_c;

d1_max = 0.05;   % hard bound on |d1| in /°C


if isMonolithic
    [log_c_ref_cold, d1_cold] = fitLogDamping(dT_c, c_v, d1_max);
    c_ref_cold = exp(log_c_ref_cold);

    c_temp     = c_ref_cold * exp(d1_cold * dT_c);
    Phi_amp    = [A_c, A_c.^2];
    b_amp      = Phi_amp \ (c_v - c_temp);
    d2_cold = b_amp(1);   d3_cold = b_amp(2);
    c_ref_warm = c_ref_cold;  d1_warm = d1_cold;
    d2_warm    = d2_cold;     d3_warm = d3_cold;

    c_pred_pw = c_ref_cold*exp(d1_cold*dT_c) + d2_cold*A_c + d3_cold*A_c.^2;
    R2_c_pw   = 1 - sum((c_v-c_pred_pw).^2)/sum((c_v-mean(c_v)).^2);
    R2_c_v3   = R2_c_pw;

    fprintf('Single-regime: c_ref=%.4f  d1=%.4e /°C  d2=%.4e  d3=%.4e\n', ...
        c_ref_cold, d1_cold, d2_cold, d3_cold);
    fprintf('R²: old=%.4f  new(exp)=%.4f\n', R2_c_old, R2_c_pw);

else  % Jointed — piecewise regimes
    cold_idx = cold_c == 1;
    warm_idx = warm_c == 1;

    [log_c_ref_cold, d1_cold] = fitLogDamping(dT_c(cold_idx), c_v(cold_idx), d1_max);
    [log_c_ref_warm, d1_warm] = fitLogDamping(dT_c(warm_idx), c_v(warm_idx), d1_max);
    c_ref_cold = exp(log_c_ref_cold);
    c_ref_warm = exp(log_c_ref_warm);

    b_c = [A_c(cold_idx), A_c(cold_idx).^2] \ (c_v(cold_idx) - c_ref_cold*exp(d1_cold*dT_c(cold_idx)));
    b_w = [A_c(warm_idx), A_c(warm_idx).^2] \ (c_v(warm_idx) - c_ref_warm*exp(d1_warm*dT_c(warm_idx)));
    d2_cold = b_c(1);  d3_cold = b_c(2);
    d2_warm = b_w(1);  d3_warm = b_w(2);

    c_temp = zeros(size(c_v));
    c_temp(cold_idx) = c_ref_cold*exp(d1_cold*dT_c(cold_idx));
    c_temp(warm_idx) = c_ref_warm*exp(d1_warm*dT_c(warm_idx));
    c_pred_pw = c_temp + cold_c.*(d2_cold*A_c + d3_cold*A_c.^2) ...
                       + warm_c.*(d2_warm*A_c + d3_warm*A_c.^2);
    R2_c_pw = 1 - sum((c_v-c_pred_pw).^2)/sum((c_v-mean(c_v)).^2);

    Phi_c_v3  = [cold_c, cold_c.*dT_c, cold_c.*A_c, warm_c, warm_c.*dT_c, warm_c.*A_c];
    c_pred_v3 = Phi_c_v3 * (Phi_c_v3 \ c_v);
    R2_c_v3   = 1 - sum((c_v-c_pred_v3).^2)/sum((c_v-mean(c_v)).^2);

    fprintf('Cold: c_ref=%.4f  d1=%.4e /°C  d2=%.4e  d3=%.4e\n', c_ref_cold,d1_cold,d2_cold,d3_cold);
    fprintf('Warm: c_ref=%.4f  d1=%.4e /°C  d2=%.4e  d3=%.4e\n', c_ref_warm,d1_warm,d2_warm,d3_warm);
    fprintf('R²: old=%.4f  v3=%.4f  new(exp)=%.4f\n', R2_c_old, R2_c_v3, R2_c_pw);
end

d1_cold_n = d1_cold;   % /°C on log scale — no further normalisation
d2_cold_n = d2_cold/(c_ref_cold+eps);
d3_cold_n = d3_cold/(c_ref_cold+eps);
d1_warm_n = d1_warm;
d2_warm_n = d2_warm/(c_ref_warm+eps);
d3_warm_n = d3_warm/(c_ref_warm+eps);

% -------------------------------------------------------------------------
% BFS-BASED DESCRIPTIVE RESPONSE CLASSIFICATION
% Uses BFS sign (positive = hardening, negative = softening) from measured
% resonance frequencies. Threshold ±0.10% separates signal from noise.
% A secondary SNR check on deltaK retains Weak/Indeterminate for genuinely
% noisy temperature levels.
% -------------------------------------------------------------------------
BFS_thresh_pct = 0.10;

for Z = 1:nT
    f_pts   = proc.steppedsine.backbone{Z}.f_pts;
    valid_f = ~isnan(f_pts);

    if sum(valid_f) < 2
        summaryNonlinearityType{Z} = 'Indeterminate';
        continue;
    end

    f_valid = f_pts(valid_f);
    bfs_pct = (f_valid(end) - f_valid(1)) / f_valid(1) * 100;

    dK_vec   = arrayfun(@(a) proc.steppedsine.amp{Z}(a).Stiffness.deltaK, 1:Na);
    dK_clean = dK_vec(~isnan(dK_vec));
    dK_snr   = 0;
    if numel(dK_clean) >= 3
        dK_snr = abs(mean(dK_clean)) / (std(dK_clean) + eps);
    end

    if abs(bfs_pct) < BFS_thresh_pct && dK_snr < 3.0
        summaryNonlinearityType{Z} = 'Weak/Indeterminate';
    elseif bfs_pct >= BFS_thresh_pct
        summaryNonlinearityType{Z} = 'Hardening';
    elseif bfs_pct <= -BFS_thresh_pct
        summaryNonlinearityType{Z} = 'Softening';
    else
        summaryNonlinearityType{Z} = 'Weak/Indeterminate';
    end
end

% Post-hoc model-based onset
summaryOnset_model   = NaN(1, nT);
summaryOnsetUm_model = NaN(1, nT);
for Z = 1:nT
    g0 = gamma0_cold;
    if T_vec(Z)-T0 >= T_star_opt; g0 = gamma0_warm; end
    prev_exceed = false;
    for a = 1:Na
        A_m = AMatrix(a,Z);
        if ~isnan(A_m) && abs(g0 * A_m^2) * 100 > onset_thresh_pct
            if prev_exceed
                summaryOnset_model(Z)   = config.steppedsine.amplitude.Vector(a-1);
                summaryOnsetUm_model(Z) = AMatrix(a-1, Z) * 1e6;
                break;
            end
            prev_exceed = true;
        else
            prev_exceed = false;
        end
    end
    if isnan(summaryOnset_model(Z))
        summaryOnset_model(Z)   = config.steppedsine.amplitude.Vector(end);
        summaryOnsetUm_model(Z) = AMatrix(Na, Z) * 1e6;
    end
end

% Store results
proc.thermo.stiffness.f_ref_Hz = f_ref_Hz;
proc.thermo.stiffness.coeffs = struct('k_ref',k_ref_fit,'c1',c1_fit, ...
    'gamma0',gamma0_fit,'gamma1',gamma1_fit, ...
    'gamma0_cold',gamma0_cold,'gamma1_cold',gamma1_cold, ...
    'gamma0_warm',gamma0_warm,'gamma1_warm',gamma1_warm, ...
    'T_star_dT',T_star_opt,'T_star_abs',T0+T_star_opt, ...
    'T_star_fixed',use_fixed_Tstar);
proc.thermo.stiffness.R2              = R2_k;
proc.thermo.stiffness.R2_global       = R2_k_global;
proc.thermo.stiffness.R2_sweep        = R2_k_sweep;
proc.thermo.stiffness.T_star_sweep    = T_star_candidates;
proc.thermo.stiffness.A_data          = A_v;
proc.thermo.stiffness.dT_data         = dT_v;
proc.thermo.stiffness.k_meas          = k_v;
proc.thermo.stiffness.k_pred          = k_pred;
proc.thermo.stiffness.resid_by_T_mean = T_resid_mean;
proc.thermo.stiffness.resid_by_T_std  = T_resid_std;
proc.thermo.stiffness.resid_by_A_mean = A_resid_mean;
proc.thermo.stiffness.resid_by_A_std  = A_resid_std;

proc.thermo.damping.coeffs = struct('c_ref_cold',c_ref_cold,'c_ref_warm',c_ref_warm, ...
    'd1_cold',d1_cold,'d1_warm',d1_warm,'d2_cold',d2_cold,'d2_warm',d2_warm, ...
    'd3_cold',d3_cold,'d3_warm',d3_warm,'d1_cold_n',d1_cold_n,'d1_warm_n',d1_warm_n, ...
    'd2_cold_n',d2_cold_n,'d2_warm_n',d2_warm_n,'d3_cold_n',d3_cold_n,'d3_warm_n',d3_warm_n, ...
    'T_star_abs',T0+T_star_opt,'MAD_thresh',MAD_thresh,'n_outliers',n_flagged, ...
    'd1_bound',d1_max,'damping_model','exp_bounded');
proc.thermo.damping.R2           = R2_c_pw;
proc.thermo.damping.R2_v3        = R2_c_v3;
proc.thermo.damping.R2_old       = R2_c_old;
proc.thermo.damping.c_meas       = c_v;
proc.thermo.damping.c_pred       = c_pred_pw;
proc.thermo.damping.regime       = regime_c;
proc.thermo.damping.outlier_flag = c_outlier_flag;
proc.thermo.damping.c_meas_all   = c_v_all;
proc.thermo.damping.A_all        = A_c_all;
proc.thermo.damping.dT_all       = dT_c_all;
proc.thermo.onset.loop           = summaryOnset;
proc.thermo.onset.loop_um        = summaryOnsetUm;
proc.thermo.onset.model          = summaryOnset_model;
proc.thermo.onset.model_um       = summaryOnsetUm_model;

% =========================================================================
% CSV EXPORT
% =========================================================================

switch config.steppedsine.system.BeamType
    case 2;    systemLabel = 'Monolithic';
    case 3;    systemLabel = 'Jointed';
    otherwise; systemLabel = sprintf('BeamType%d', config.steppedsine.system.BeamType);
end

outlier_lookup = false(Na, nT);
for ai_ol = 1:Na
    for zi_ol = 1:nT
        flat_ol = (zi_ol-1)*Na + ai_ol;
        if flat_ol <= numel(valid_c_raw) && valid_c_raw(flat_ol)
            outlier_lookup(ai_ol,zi_ol) = c_outlier_flag(sum(valid_c_raw(1:flat_ol)));
        end
    end
end

exportRows = {};
for zi = 1:nT
    T_val  = config.steppedsine.temperature.Vector(zi);
    dT_val = T_val - T_vec(1);
    for ai = 1:Na
        exportRows{end+1} = {systemLabel, T_val, dT_val, ...
            config.steppedsine.amplitude.Vector(ai), AMatrix(ai,zi), kEffMatrix(ai,zi), ...
            (kEffMatrix(ai,zi)-k_ref_fit)/k_ref_fit, zetaMatrix(ai,zi), cEqMatrix(ai,zi), ...
            peakFreqMatrix(ai,zi), peakAmpMatrix(ai,zi), ...
            double(outlier_lookup(ai,zi)), proc.steppedsine.amp{zi}(ai).DampingMethod};
    end
end

obsTable = cell2table(vertcat(exportRows{:}), 'VariableNames', ...
    {'SystemType','Temperature_degC','DeltaT_degC','ExcitationAmplitude', ...
     'A_metres','k_eff_Npm','deltaK_normalised','Zeta','c_eq_Nspm', ...
     'ResonanceFreq_Hz','PeakAcc_g','DampingOutlierFlag','DampingMethod'});

% Export fields for the retained piecewise thermal model
coeffTable = table({systemLabel}, k_ref_fit, f_ref_Hz, c1_fit, ...
    R2_k, R2_k_global, T0+T_star_opt, double(use_fixed_Tstar), ...
    gamma0_cold, gamma1_cold, gamma0_warm, gamma1_warm, ...
    c_ref_cold, c_ref_warm, d1_cold, d1_warm, d2_cold, d2_warm, d3_cold, d3_warm, ...
    n_flagged, R2_c_pw, R2_c_v3, R2_c_old, ...
    'VariableNames', {'SystemType','k_ref_Npm','f_ref_Hz','c1_per_degC', ...
    'R2_stiffness_piecewise','R2_stiffness_global','T_star_degC','T_star_fixed', ...
    'gamma0_cold_per_m2','gamma1_cold_per_m2_degC','gamma0_warm_per_m2','gamma1_warm_per_m2_degC', ...
    'c_ref_cold_Nspm','c_ref_warm_Nspm','d1_cold_per_degC','d1_warm_per_degC', ...
    'd2_cold_Nspm_per_m','d2_warm_Nspm_per_m','d3_cold_Nspm_per_m2','d3_warm_Nspm_per_m2', ...
    'n_damping_outliers','R2_damping','R2_damping_v3','R2_damping_old'});

summaryRows = {};
for zi = 1:nT
    T_val = config.steppedsine.temperature.Vector(zi);
    bb    = proc.steppedsine.backbone{zi};
    fShift = NaN;
    if ~all(isnan(bb.f_pts))
        fpts = bb.f_pts(~isnan(bb.f_pts));
        if numel(fpts) >= 2; fShift = (fpts(end)-fpts(1))/fpts(1); end
    end
    zetas = arrayfun(@(a) proc.steppedsine.amp{zi}(a).Zeta, 1:Na);
    summaryRows{end+1} = {systemLabel, T_val, summaryOnset(zi), summaryOnsetUm(zi), ...
                           summaryNonlinearityType{zi}, fShift, mean(zetas,'omitnan')};
end
summaryTable = cell2table(vertcat(summaryRows{:}), 'VariableNames', ...
    {'SystemType','Temperature_degC','OnsetAmplitude_V','OnsetDisplacement_um', ...
     'NonlinearityType','BackboneFreqShift','MeanZeta'});

if ~exist('exportCache','var'); exportCache = struct(); end
exportCache.(systemLabel).obs     = obsTable;
exportCache.(systemLabel).coeffs  = coeffTable;
exportCache.(systemLabel).summary = summaryTable;

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
writetable(obsTable,     sprintf('export_%s_observations_%s.csv',  systemLabel, timestamp));
writetable(coeffTable,   sprintf('export_%s_coefficients_%s.csv',  systemLabel, timestamp));
writetable(summaryTable, sprintf('export_%s_summary_%s.csv',       systemLabel, timestamp));


% Power law export (jointed only — NaN rows included for completeness)
T_export = config.steppedsine.temperature.Vector(:);
powerlaw_table = table(T_export, powerlaw_n(:), powerlaw_C(:), powerlaw_R2(:), ...
    'VariableNames', {'Temperature_degC', 'Mindlin_n', 'Mindlin_C', 'Mindlin_R2'});
writetable(powerlaw_table, ...
    sprintf('export_%s_powerlaw_%s.csv', systemLabel, timestamp));
fprintf('Exported: power law n values — %d rows\n', height(powerlaw_table));

fprintf('Exported: %s obs=%d rows\n', systemLabel, height(obsTable));




% % =========================================================================
% % INTERACTIVE TEMPERATURE BROWSER
% % =========================================================================
%
% switch config.steppedsine.system.BeamType
%     case 2; figName = 'Temperature Browser — Monolithic';
%     case 3; figName = 'Temperature Browser — Jointed';
%     otherwise; figName = 'Temperature Browser';
% end
% fig = figure('Name', figName);
% ax  = axes(fig); hold(ax,'on'); grid(ax,'on');
% xlabel(ax,'Frequency (Hz)'); ylabel(ax,'Response Amplitude');
% set(ax,'YScale','log');
% xlim(sort([config.steppedsine.freq.Start, config.steppedsine.freq.End]));
% ylim([0 100]);
%
% setappdata(fig,'activeZ', 1);
% setappdata(fig,'temps',   config.steppedsine.temperature.Vector);
% setappdata(fig,'Na',      Na);
% setappdata(fig,'proc',    proc);
% setappdata(fig,'ax',      ax);
%
% drawTemperature(fig);
% set(fig,'KeyPressFcn',@temperatureCycleCallback);

% ╔══════════════════════════════════════════════════════════════════════╗
% ║ INTERACTIVE TEMPERATURE BROWSER                                     ║
% ╚══════════════════════════════════════════════════════════════════════╝

if config.steppedsine.system.BeamType == 2
    fig = figure('Name','Temperature Browser - Monolithic');
elseif config.steppedsine.system.BeamType == 3
    fig = figure('Name','Temperature Browser - Jointed');
else
    error('ERROR')
end
ax  = axes(fig);
hold(ax,'on'); grid(ax,'on');

xlabel(ax,'Frequency (Hz)');
ylabel(ax,'Response Amplitude');
set(ax,'YScale','log');
xlim([config.steppedsine.freq.Start config.steppedsine.freq.End]); ylim([0 100]);

nT    = numel(config.steppedsine.temperature.Vector);

% Storage for plot handles
hAmp  = gobjects(Na,1);
hBB   = gobjects(1);

activeZ = 1;

% Store shared data for callback
setappdata(fig,'activeZ',activeZ);
setappdata(fig,'temps',config.steppedsine.temperature.Vector);
setappdata(fig,'Na',Na);
setappdata(fig,'proc',proc);
setappdata(fig,'ax',ax);

% ============================================================
% SECOND FIGURE — FRF BROWSER (independent window)
% ============================================================

figFRF = figure('Name','FRF Browser');
axFRF  = axes(figFRF);
hold(axFRF,'on'); grid(axFRF,'on');

xlabel(axFRF,'Frequency (Hz)');
ylabel(axFRF,'|H(f)|');
set(axFRF,'YScale','log');

xlim(axFRF,[config.steppedsine.freq.Start config.steppedsine.freq.End]);

% Store required data
setappdata(figFRF,'proc',proc);
setappdata(figFRF,'temps',config.steppedsine.temperature.Vector);
setappdata(figFRF,'Na',Na);
setappdata(figFRF,'ax',axFRF);

% Link both figures together
setappdata(fig,'figFRF',figFRF);
setappdata(figFRF,'figMain',fig);

% Initial draw
drawTemperature(fig);

% --- Update FRF browser if it exists ---
if isappdata(fig,'figFRF')
    figFRF = getappdata(fig,'figFRF');
elseif isappdata(fig,'figMain')
    figFRF = fig;
    fig    = getappdata(fig,'figMain');
else
    return
end

if isvalid(figFRF)
    drawFRF(figFRF);
end

% Key control
set(fig,'KeyPressFcn',@temperatureCycleCallback);
set(figFRF,'KeyPressFcn',@temperatureCycleCallback);

% =========================================================================
% INDIVIDUAL BACKBONE PLOTS — one per temperature, using browser's draw logic
% =========================================================================

bbDir = fullfile('exports', 'figures', 'backbones', systemLabel);
if ~exist(bbDir, 'dir'); mkdir(bbDir); end

for zi = 1:nT

    T_val = config.steppedsine.temperature.Vector(zi);
    bb    = proc.steppedsine.backbone{zi};

    figBB = figure('Name', sprintf('%s  T = %d°C', systemLabel, T_val), ...
                   'Visible', 'off');
    ax_bb = axes(figBB);
    hold(ax_bb, 'on');
    grid(ax_bb, 'on');
    set(ax_bb, 'YScale', 'log');
    xlim(ax_bb, sort([config.steppedsine.freq.Start, config.steppedsine.freq.End]));
    ylim(ax_bb, [0 100]);
    xlabel(ax_bb, 'Frequency (Hz)');
    ylabel(ax_bb, 'Response Amplitude');
    % title(ax_bb, sprintf('%s — T = %d°C', systemLabel, T_val));

    % Replicate drawTemperature: plot each amplitude's peak-fit curve
    for a = 1:Na
        amp_a = proc.steppedsine.amp{zi}(a);
        if ~isempty(amp_a.Fit.fDense) && ~isempty(amp_a.Fit.yFit)
            plot(ax_bb, amp_a.Fit.fDense, amp_a.Fit.yFit, 'b-', 'LineWidth', 1);
        end
    end

    % % Backbone points
    % validBB = ~isnan(bb.f_pts) & ~isnan(bb.y_pts);
    % if any(validBB)
    %     scatter(ax_bb, bb.f_pts(validBB), bb.y_pts(validBB), 40, 'r', ...
    %             'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    % end

    % Backbone curve
    if ~all(isnan(bb.f)) && ~all(isnan(bb.y))
        plot(ax_bb, bb.f, bb.y, 'r-', 'LineWidth', 2);
    end

    exportgraphics(figBB, ...
        fullfile(bbDir, sprintf('backbone_%s_T%+04d.png', systemLabel, T_val)), ...
        'Resolution', 150);
    % close(figBB);

end
fprintf('Backbone figures saved to %s\n', bbDir);

fprintf('Backbone figures saved to %s\n', bbDir);

if config.steppedsine.system.BeamType == 2
    proc.backbone.monolithic = proc.steppedsine.backbone;
else
    proc.backbone.jointed = proc.steppedsine.backbone;
end
%% ==== END SteppedSine2.m ====

proc_all.Jointed.steppedsine = proc.steppedsine;
proc_all.Jointed.thermo      = proc.thermo;
% hysteresis.Jointed.upsweep	= proc.steppedsine;
% config.steppedsine.data.fileName = arrayfun(@(k) ...
%     sprintf("0001-00%02d_117-112-0.1_1000_2000_3.csv", k), 1:config.steppedsine.nTemps);
% run SteppedSine2.m
% hysteresis.Jointed.downsweep	= proc.steppedsine;

% Combined CSV export (done here, once, after both systems are processed)
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
if isfield(proc_all.Monolithic,'thermo') && isfield(proc_all.Jointed,'thermo')
    % Pull the per-system obs/coeffs/summary tables saved by SteppedSine2
    % (they're in exportCache which SteppedSine2 populates in the caller workspace)
    if exist('exportCache','var') && isfield(exportCache,'Monolithic') && isfield(exportCache,'Jointed')
        combined_obs     = [exportCache.Monolithic.obs;     exportCache.Jointed.obs];
        combined_coeffs  = [exportCache.Monolithic.coeffs;  exportCache.Jointed.coeffs];
        combined_summary = [exportCache.Monolithic.summary; exportCache.Jointed.summary];
        writetable(combined_obs,     sprintf('export_COMBINED_observations_%s.csv',  timestamp));
        writetable(combined_coeffs,  sprintf('export_COMBINED_coefficients_%s.csv',  timestamp));
        writetable(combined_summary, sprintf('export_COMBINED_summary_%s.csv',       timestamp));

        mono_t = exportCache.Monolithic.obs;
        jnt_t  = exportCache.Jointed.obs;
        [~,iM,iJ] = intersect([mono_t.Temperature_degC, mono_t.ExcitationAmplitude], ...
                               [jnt_t.Temperature_degC,  jnt_t.ExcitationAmplitude], 'rows');
        if ~isempty(iM)
            compTable = table(mono_t.Temperature_degC(iM), mono_t.ExcitationAmplitude(iM), ...
                mono_t.k_eff_Npm(iM), jnt_t.k_eff_Npm(iJ), ...
                jnt_t.k_eff_Npm(iJ)./mono_t.k_eff_Npm(iM), ...
                jnt_t.k_eff_Npm(iJ)-mono_t.k_eff_Npm(iM), ...
                jnt_t.Zeta(iJ)-mono_t.Zeta(iM), ...
                jnt_t.deltaK_normalised(iJ)-mono_t.deltaK_normalised(iM), ...
                'VariableNames',{'Temperature_degC','ExcitationAmplitude', ...
                'k_Monolithic_Npm','k_Jointed_Npm','k_Ratio','k_Loss_Npm', ...
                'DeltaZeta','DeltaDeltaK'});
            writetable(compTable, sprintf('export_COMPARISON_joint_effect_%s.csv', timestamp));
            fprintf('Combined export written (%d matched pairs).\n', height(compTable));
        end
    end
end

% ── Step 3: System comparison ─────────────────────────────────────────────
fprintf('\nStep 3: Comparison\n');

%% ==== BEGIN Comparison.m (integrated) ====
% Comparison.m
% Joint effect characterisation from stepped-sine backbones.
% Uses displacement amplitude (metres), normalised Area metric,
% and raw-point curvature for acrit.

fM0 = proc.impact.monolithic.peakFreq;
fJ0 = proc.impact.jointed.peakFreq;

nTemps = 16;
temps  = -30:10:120;

JK_all       = cell(nTemps,1);
acrit_all    = NaN(nTemps,1);
acrit_method = cell(nTemps,1);
Area_all     = NaN(nTemps,1);
Area_density = NaN(nTemps,1);
Awin_lo      = NaN(nTemps,1);
Awin_hi      = NaN(nTemps,1);
Astar_all    = cell(nTemps,1);
Kratio_all   = cell(nTemps,1);
JKstar_all   = cell(nTemps,1);
DeltaF_all   = cell(nTemps,1);

for i = 1:nTemps

    bbJ = proc.backbone.jointed{i};
    bbM = proc.backbone.monolithic{i};

    % Use physical displacement amplitude if available (SteppedSine2 v5+),
    % otherwise fall back to polynomial backbone acceleration axis.
    has_Apts = isfield(bbJ,'A_pts') && isfield(bbM,'A_pts') && ...
               ~all(isnan(bbJ.A_pts)) && ~all(isnan(bbM.A_pts));

    if has_Apts
        validJ = ~isnan(bbJ.A_pts) & ~isnan(bbJ.f_pts) & bbJ.A_pts>0;
        validM = ~isnan(bbM.A_pts) & ~isnan(bbM.f_pts) & bbM.A_pts>0;
        [aJ_raw,sJ] = sort(bbJ.A_pts(validJ));  fJ_raw = bbJ.f_pts(validJ); fJ_raw = fJ_raw(sJ)/fJ0;
        [aM_raw,sM] = sort(bbM.A_pts(validM));  fM_raw = bbM.f_pts(validM); fM_raw = fM_raw(sM)/fM0;
        amp_label = 'Displacement amplitude (m)';
    else
        warning('Comparison: A_pts missing for T index %d — using polynomial axis.', i);
        aJ_raw = bbJ.y(:);  fJ_raw = bbJ.f(:)/fJ0;
        aM_raw = bbM.y(:);  fM_raw = bbM.f(:)/fM0;
        amp_label = 'Response amplitude (acc)';
    end

    a_lo = max(min(aJ_raw), min(aM_raw));
    a_hi = min(max(aJ_raw), max(aM_raw));
    if a_lo >= a_hi
        warning('Comparison: no overlapping amplitude range for T index %d.', i);
        continue;
    end

    a_common   = linspace(a_lo, a_hi, 300);
    Awin_lo(i) = a_lo;
    Awin_hi(i) = a_hi;

    fJ_i = interp1(aJ_raw, fJ_raw, a_common, 'pchip', 'extrap');
    fM_i = interp1(aM_raw, fM_raw, a_common, 'pchip', 'extrap');

    DeltaF_all{i} = struct('a',a_common, ...
        'J',(fJ_i-fJ_i(1))/fJ_i(1), 'M',(fM_i-fM_i(1))/fM_i(1), 'amp_label',amp_label);

    k_ratio   = (fJ_i./fM_i).^2;
    JK        = abs(fJ_i-fM_i)./fM_i;
    JK_all{i} = JK;
    dk_da     = gradient(k_ratio, a_common);

    % Area metrics
    Area_all(i)    = trapz(a_common, abs(fJ_i-fM_i));
    Area_density(i) = Area_all(i) / (a_hi-a_lo);

    % acrit from raw discrete points (3-point finite difference)
    if has_Apts && sum(validJ)>=4
        in_win = aJ_raw>=a_lo & aJ_raw<=a_hi;
        aJ_w   = aJ_raw(in_win);
        fJ_w   = fJ_raw(in_win);
        if numel(aJ_w) >= 4
            kappa = NaN(numel(aJ_w),1);
            for ip = 2:numel(aJ_w)-1
                h_avg = (aJ_w(ip)-aJ_w(ip-1) + aJ_w(ip+1)-aJ_w(ip)) / 2;
                if h_avg>0
                    kappa(ip) = abs(fJ_w(ip+1)-2*fJ_w(ip)+fJ_w(ip-1)) / h_avg^2;
                end
            end
            [~,ik] = max(kappa,[],'omitnan');
            acrit_all(i)    = aJ_w(ik);
            acrit_method{i} = 'RawPoints';
        else
            acrit_all(i)    = a_common(find(abs(gradient(gradient(fJ_i,a_common),a_common))==max(abs(gradient(gradient(fJ_i,a_common),a_common))),1));
            acrit_method{i} = 'Polynomial';
        end
    else
        d2f = gradient(gradient(fJ_i,a_common),a_common);
        [~,idx] = max(abs(d2f));
        acrit_all(i)    = a_common(idx);
        acrit_method{i} = 'Polynomial';
    end

    a_star        = a_common / acrit_all(i);
    Astar_all{i}  = a_star;
    Kratio_all{i} = k_ratio;
    JKstar_all{i} = JK;

end

% Console summary
fprintf('\n=== JOINT EFFECT SUMMARY ===\n');
fprintf('%-6s  %-10s  %-10s  %-12s  %-12s  %-10s  %s\n', ...
    'T', 'a_lo (m)', 'a_hi (m)', 'Area', 'Area/width', 'acrit (m)', 'method');
for i = 1:nTemps
    if ~isnan(Area_all(i))
        fprintf('%-6.0f  %-10.3e  %-10.3e  %-12.4e  %-12.4e  %-10.3e  %s\n', ...
            temps(i), Awin_lo(i), Awin_hi(i), Area_all(i), Area_density(i), ...
            acrit_all(i), acrit_method{i});
    end
end
fprintf('\nArea_density (Area/window_width) is the comparable metric across temperatures.\n\n');

% ── Interactive Δf/f₀ browser ────────────────────────────────────────────
fig = figure('Name','Normalised Frequency Shift (Temperature Browser)');
ax  = axes(fig); hold(ax,'on'); grid(ax,'on');
first_valid = find(~cellfun(@isempty,DeltaF_all),1);
xlabel(ax, DeltaF_all{first_valid}.amp_label);
ylabel(ax, '\Delta f / f_0');
xlim(ax, [0, max(Awin_hi(~isnan(Awin_hi)))*1.05]);
setappdata(fig,'activeIdx',1); setappdata(fig,'temps',temps);
setappdata(fig,'DeltaF',DeltaF_all); setappdata(fig,'ax',ax);
drawDeltaF(fig);
set(fig,'KeyPressFcn',@cycleDeltaFCallback);
legend('Jointed','Monolithic');

% ── Summary figures ────────────────────────────────────────────────────────
valid_i = ~isnan(Area_all);

% figure('Name','Joint Nonlinearity Severity','Color','w');
% subplot(2,2,1);
% bar(temps(valid_i), Area_all(valid_i),'FaceColor',[0.4 0.6 0.8]);
% xlabel('T (°C)'); ylabel('Area'); title('Raw Area (window-dependent)'); grid on;
%
% subplot(2,2,2);
% bar(temps(valid_i), Area_density(valid_i),'FaceColor',[0.2 0.6 0.4]);
% xlabel('T (°C)'); ylabel('Area/width'); title('Area density (comparable across T)'); grid on;
%
% subplot(2,2,3);
% plot(temps(valid_i), Awin_lo(valid_i)*1e6,'bs-','LineWidth',1.2,'DisplayName','Window lo'); hold on;
% plot(temps(valid_i), Awin_hi(valid_i)*1e6,'rs-','LineWidth',1.2,'DisplayName','Window hi');
% xlabel('T (°C)'); ylabel('Amplitude (µm)'); title('Integration window bounds'); legend; grid on;
%
% subplot(2,2,4);
% bar(temps(valid_i), acrit_all(valid_i)*1e6,'FaceColor',[0.8 0.5 0.2]); hold on;
% fb_temps = temps(valid_i);
% fb_vals  = acrit_all(valid_i & strcmp(acrit_method,'Polynomial'))*1e6;
% if ~isempty(fb_vals)
%     plot(fb_temps(strcmp(acrit_method(valid_i),'Polynomial')), fb_vals, ...
%         'rx','MarkerSize',10,'LineWidth',2,'DisplayName','Fallback');
%     legend('a_{crit}','Polynomial fallback','Location','best');
% end
% xlabel('T (°C)'); ylabel('a_{crit} (µm)'); title('Slip onset amplitude'); grid on;
%
% % k-ratio collapse
% figure('Name','k-ratio vs Normalised Amplitude','Color','w');
% ax2=axes; hold(ax2,'on'); grid(ax2,'on');
% try; ax2.ColorOrder=tempmap; catch; end
% for i=1:nTemps
%     if ~isempty(Astar_all{i}); plot(ax2,Astar_all{i},Kratio_all{i},'LineWidth',1.2); end
% end
% xline(ax2,1,'k--','a_{crit}','LineWidth',1.2);
% xlabel(ax2,'a / a_{crit}'); ylabel(ax2,'k_{eff}/k_{mono}'); title(ax2,'k-ratio collapse'); xlim(ax2,[0 3]);
%
% % JK collapse
% figure('Name','J_K vs Normalised Amplitude','Color','w');
% ax3=axes; hold(ax3,'on'); grid(ax3,'on');
% try; ax3.ColorOrder=tempmap; catch; end
% for i=1:nTemps
%     if ~isempty(Astar_all{i}); plot(ax3,Astar_all{i},JKstar_all{i},'LineWidth',1.2); end
% end
% xline(ax3,1,'k--','a_{crit}','LineWidth',1.2);
% xlabel(ax3,'a / a_{crit}'); ylabel(ax3,'J_K'); title(ax3,'JNI metric J_K collapse'); xlim(ax3,[0 3]);
%% ==== END Comparison.m ====

% ── Step 4: Thermo-mechanical coupling ────────────────────────────────────
fprintf('\nStep 4: TMCoupling\n');

%% ==== BEGIN TMCoupling.m (integrated) ====

T0      = 20;
outDir  = fullfile('exports','figures');
if ~exist(outDir,'dir'); mkdir(outDir); end
systems = {'Monolithic','Jointed'};

if ~exist('proc_all','var'); error('proc_all not found. Run process2.m first.'); end
for sIdx = 1:2
    if ~isfield(proc_all, systems{sIdx}) || ~isfield(proc_all.(systems{sIdx}),'thermo')
        error('proc_all.%s.thermo not found. Check process2.m saves proc.thermo.', systems{sIdx});
    end
end

results = struct();

for sIdx = 1:2
    sys = systems{sIdx};
    th  = proc_all.(sys).thermo;

    A_v       = th.stiffness.A_data;
    dT_v      = th.stiffness.dT_data;
    k_v       = th.stiffness.k_meas;
    k_pred_pw = th.stiffness.k_pred;
    T_v       = dT_v + T0;
    k_ref_fit = th.stiffness.coeffs.k_ref;
    deltaK_v  = (k_v - k_ref_fit) / k_ref_fit;

    fprintf('\n=== %s ===\n', sys);
    fprintf('Piecewise: k_ref=%.2f N/m  c1=%.4e /C  T*=%.0fC\n', ...
        k_ref_fit, th.stiffness.coeffs.c1, th.stiffness.coeffs.T_star_abs);
    fprintf('Cold: g0=%.4e  g1=%.4e  /m2\n', th.stiffness.coeffs.gamma0_cold, th.stiffness.coeffs.gamma1_cold);
    fprintf('Warm: g0=%.4e  g1=%.4e  /m2\n', th.stiffness.coeffs.gamma0_warm, th.stiffness.coeffs.gamma1_warm);
    fprintf('R2: piecewise=%.4f  global=%.4f\n', th.stiffness.R2, th.stiffness.R2_global);

    % Additional polynomial diagnostic on the same data. Scale predictors to
    % O(1), then remove terms that the observations cannot distinguish.
    % (In particular, A^4 in metres^4 is ~1e-24 for micrometre motion.)
    ampScale = max(abs(A_v(isfinite(A_v))));
    tempScale = max(abs(dT_v(isfinite(dT_v))));
    if isempty(ampScale) || ampScale == 0, ampScale = 1; end
    if isempty(tempScale) || tempScale == 0, tempScale = 1; end
    tbl = table(dT_v(:)/tempScale, A_v(:)/ampScale, deltaK_v(:), ...
        'VariableNames', {'dT','A','deltaK'});
    tbl.A2   = tbl.A.^2;   tbl.A4   = tbl.A.^4;
    tbl.dTA2 = tbl.dT.*tbl.A2;   tbl.dTA4 = tbl.dT.*tbl.A4;
    terms = {'dT','A2','A4','dTA2','dTA4'};
    usable = all(isfinite(tbl{:, [{'deltaK'}, terms]}), 2);
    if nnz(usable) < 3
        error('Additional stiffness-related regression requires at least 3 finite observations for %s.', sys);
    end
    % Greedy rank test includes the intercept and preserves the scientific
    % hierarchy (linear terms before higher-order and interaction terms).
    X = ones(nnz(usable),1);
    selectedTerms = {};
    for it = 1:numel(terms)
        candidate = [X, tbl{usable,terms{it}}];
        if rank(candidate, 1e-9) > size(X,2)
            X = candidate;
            selectedTerms{end+1} = terms{it}; %#ok<SAGROW>
        else
            fprintf('  Unidentifiable additional term omitted: %s (%s)\n', terms{it}, sys);
        end
    end
    if isempty(selectedTerms)
        modelFormula = 'deltaK ~ 1';
    else
        modelFormula = ['deltaK ~ ', strjoin(selectedTerms, ' + ')];
    end
    mdl = fitlm(tbl(usable,:), modelFormula);
    fprintf('  Regression coefficients below use normalised temperature/amplitude predictors.\n');
    fprintf('Polynomial (additional diagnostic): R2=%.4f  Adj R2=%.4f\n', ...
        mdl.Rsquared.Ordinary, mdl.Rsquared.Adjusted);

    ci = coefCI(mdl, 0.05);
    fprintf('  %-20s  %+12s  %-20s  %-8s\n','Coefficient','Estimate','95% CI','p-value');
    for ii = 1:numel(mdl.CoefficientNames)
        p = mdl.Coefficients.pValue(ii);
        if p<0.001; sig='***'; elseif p<0.01; sig='**'; elseif p<0.05; sig='*'; else; sig='(ns)'; end
        fprintf('  %-20s  %+12.4e  [%+.3e %+.3e]  %.4f %s\n', ...
            mdl.CoefficientNames{ii}, mdl.Coefficients.Estimate(ii), ci(ii,1), ci(ii,2), p, sig);
    end

    % Restore coefficients to the original (Celsius, metres) variables
    % before predicting surfaces or comparing with piecewise coefficients.
    C_int = getVar(mdl,'(Intercept)');
    C1    = getVar(mdl,'dT')   / tempScale;
    C3    = getVar(mdl,'A2')   / ampScale^2;
    C4    = getVar(mdl,'A4')   / ampScale^4;
    C6    = getVar(mdl,'dTA2') / (tempScale*ampScale^2);
    C7    = getVar(mdl,'dTA4') / (tempScale*ampScale^4);

    if abs(C1 - th.stiffness.coeffs.c1) / (abs(th.stiffness.coeffs.c1)+eps) > 0.25
        fprintf('  WARNING: polynomial C1=%.4e vs piecewise c1=%.4e differ >25%%\n', C1, th.stiffness.coeffs.c1);
    end

    % Per-temperature A^2 slope
    uniqueT = unique(T_v);
    eta1    = NaN(size(uniqueT));
    for iT = 1:numel(uniqueT)
        idx_T = T_v == uniqueT(iT);
        if sum(idx_T) >= 3
            p_eta = polyfit(A_v(idx_T).^2, deltaK_v(idx_T), 1);
            eta1(iT) = p_eta(1);
        end
    end

    % DeltaK surfaces for plotting
    [Agrid,Tgrid] = meshgrid(linspace(min(A_v),max(A_v),50), linspace(min(T_v),max(T_v),50));
    dTg = Tgrid - T0;
    r_g = double(dTg >= th.stiffness.coeffs.T_star_abs - T0);

    deltaKg_poly = C_int + C1*dTg + C3*Agrid.^2 + C4*Agrid.^4 + C6*dTg.*Agrid.^2 + C7*dTg.*Agrid.^4;
    deltaKg_pw   = th.stiffness.coeffs.c1*dTg + ...
        (1-r_g).*(th.stiffness.coeffs.gamma0_cold*Agrid.^2 + th.stiffness.coeffs.gamma1_cold*dTg.*Agrid.^2) + ...
        r_g.*    (th.stiffness.coeffs.gamma0_warm*Agrid.^2 + th.stiffness.coeffs.gamma1_warm*dTg.*Agrid.^2);

    % Store
    results.(sys) = struct('k_ref',k_ref_fit,'c1_pw',th.stiffness.coeffs.c1,'c1_poly',C1, ...
        'C3',C3,'C4',C4,'C6',C6,'R2_pw',th.stiffness.R2, ...
        'R2_poly',mdl.Rsquared.Ordinary,'R2_adj_poly',mdl.Rsquared.Adjusted, ...
        'uniqueT',uniqueT,'eta1',eta1,'A_v',A_v,'T_v',T_v,'deltaK_v',deltaK_v, ...
        'deltaK_pred',(k_pred_pw-k_ref_fit)/k_ref_fit,'mdl',mdl, ...
        'Agrid',Agrid,'Tgrid',Tgrid,'deltaKg_poly',deltaKg_poly,'deltaKg_pw',deltaKg_pw);

    % Per-system figures
%     figure('Color','w','Name',sprintf('%s — DeltaK surfaces',sys));
%     subplot(1,2,1);
%     surf(Agrid*1e6,Tgrid,deltaKg_poly,'FaceAlpha',0.75,'EdgeColor','none'); hold on;
%     scatter3(A_v*1e6,T_v,deltaK_v,24,'k','filled');
%     xlabel('A (µm)'); ylabel('T (°C)'); zlabel('\DeltaK'); title('Polynomial'); grid on; view(35,30);
%     subplot(1,2,2);
%     surf(Agrid*1e6,Tgrid,deltaKg_pw,'FaceAlpha',0.75,'EdgeColor','none'); hold on;
%     scatter3(A_v*1e6,T_v,deltaK_v,24,'k','filled');
%     xlabel('A (µm)'); ylabel('T (°C)'); zlabel('\DeltaK');
%     title(sprintf('Piecewise T*=%.0f°C',th.stiffness.coeffs.T_star_abs)); grid on; view(35,30);
%     sgtitle(sprintf('%s: \\DeltaK surfaces',sys));
%     saveas(gcf, fullfile(outDir,sprintf('%s_DeltaK_surfaces.png',sys)));
%
%     figure('Color','w','Name',sprintf('%s — Residuals',sys));
%     deltaK_pred_v = (k_pred_pw - k_ref_fit)/k_ref_fit;
%     res_poly = mdl.Residuals.Raw;
%     subplot(2,2,1); scatter(deltaK_pred_v, deltaK_v-deltaK_pred_v, 15, T_v,'filled');
%     colorbar; yline(0,'k--'); xlabel('Predicted \DeltaK'); ylabel('Residual'); title('Piecewise'); grid on;
%     subplot(2,2,2); scatter(deltaK_v, res_poly, 15, T_v,'filled');
%     colorbar; yline(0,'k--'); xlabel('Measured \DeltaK'); ylabel('Poly residual'); title('Polynomial'); grid on;
%     subplot(2,2,3); qqplot(deltaK_v-deltaK_pred_v); title('Piecewise normality');
%     subplot(2,2,4); qqplot(res_poly); title('Polynomial normality');
%     sgtitle(sprintf('%s: Residuals',sys));
%     saveas(gcf, fullfile(outDir,sprintf('%s_Residuals.png',sys)));
%
%     figure('Color','w','Name',sprintf('%s — eta1 vs T',sys));
%     plot(uniqueT, eta1,'o-','LineWidth',1.5,'MarkerFaceColor','auto');
%     xline(th.stiffness.coeffs.T_star_abs,'r--','Label',sprintf('T*=%.0f°C',th.stiffness.coeffs.T_star_abs));
%     xlabel('T (°C)'); ylabel('\eta_1  (1/m^2)'); title(sprintf('%s: A^2 slope vs T',sys)); grid on;
%     saveas(gcf, fullfile(outDir,sprintf('%s_eta1_vs_T.png',sys)));
end

% Comparison table
fprintf('\n=== COMPARISON ===\n');
fprintf('%-28s  %-14s  %-14s\n','Parameter','Monolithic','Jointed');
fprintf('%s\n', repmat('-',1,58));
rows = {'k_ref (N/m)','k_ref','%.2f'; 'c1 piecewise (/C)','c1_pw','%.4e'; ...
        'c1 polynomial (/C)','c1_poly','%.4e'; 'C3 (A^2)','C3','%.4e'; ...
        'C6 (dT*A^2)','C6','%.4e'; 'R2 piecewise','R2_pw','%.4f'; ...
        'R2 polynomial','R2_poly','%.4f'; 'R2 adj polynomial','R2_adj_poly','%.4f'};
for ii = 1:size(rows,1)
    vM = results.Monolithic.(rows{ii,2});
    vJ = results.Jointed.(rows{ii,2});
    fmt = sprintf('%%- 28s  %%%s  %%%s\\n', rows{ii,3}(2:end), rows{ii,3}(2:end));
    fprintf(fmt, rows{ii,1}, vM, vJ);
end

% Cross-system stiffness comparison using resonance frequency ratio
% (mass-independent; k proportional to f^2 when masses are equal)
fprintf('\n=== CROSS-SYSTEM STIFFNESS COMPARISON (frequency-based) ===\n');
if isfield(proc_all.Monolithic.thermo.stiffness,'f_ref_Hz') && ...
   isfield(proc_all.Jointed.thermo.stiffness,'f_ref_Hz')
    fM = proc_all.Monolithic.thermo.stiffness.f_ref_Hz;
    fJ = proc_all.Jointed.thermo.stiffness.f_ref_Hz;
    k_ratio_freq = (fJ/fM)^2;
    fprintf('  f_Monolithic = %.4f Hz\n', fM);
    fprintf('  f_Jointed    = %.4f Hz\n', fJ);
    fprintf('  (f_J/f_M)^2  = %.4f\n', k_ratio_freq);
    fprintf('  Joint stiffness = %.1f%% of monolithic (joint reduces stiffness by %.1f%%)\n', ...
        100*k_ratio_freq, 100*(1-k_ratio_freq));
    fprintf('  NOTE: k_ref values are physical-mass-based stiffness-related scales,\n');
    fprintf('        not directly comparable across systems. Use f^2 ratio above.\n\n');
else
    fprintf('  f_ref_Hz not found in proc.thermo — re-run SteppedSine2 to populate.\n\n');
end

% Overlay eta1 comparison
figure('Color','w','Name','eta1 comparison — both systems');
hold on; grid on;
plot(results.Monolithic.uniqueT, results.Monolithic.eta1,'bs-','LineWidth',1.5,'MarkerFaceColor','b','DisplayName','Monolithic');
plot(results.Jointed.uniqueT,    results.Jointed.eta1,   'rs-','LineWidth',1.5,'MarkerFaceColor','r','DisplayName','Jointed');
for sIdx = 1:2
    T_s = proc_all.(systems{sIdx}).thermo.stiffness.coeffs.T_star_abs;
    col = 'b'; if sIdx==2; col='r'; end
    xline(T_s,[col '--'],'Label',sprintf('T*_{%s}=%.0f°C',systems{sIdx},T_s),'LabelHorizontalAlignment','left');
end
xlabel('T (°C)'); ylabel('\eta_1  (1/m^2)'); title('A^2 nonlinearity slope vs T'); legend;
saveas(gcf, fullfile(outDir,'COMPARISON_eta1_vs_T.png'));

fprintf('\nTMCoupling complete.\n\n');

% =========================================================================
% LOCAL FUNCTION
% =========================================================================
%% ==== END TMCoupling.m ====

% % ── Step 5: Hysteresis analysis ───────────────────────────────────────────
% fprintf('\nStep 5: Hysteresis\n');
% run Hysteresis.m
%
% % ── Step 6: Hysteresis plots ──────────────────────────────────────────────
% fprintf('\nStep 6: Hysteresis plots\n');
% config.hysteresis.csv = sprintf('export_Hysteresis_%s.csv', timestamp);
% run Plothysteresis.m

% ── Step 7: Duffing coefficient plot ──────────────────────────────────────
fprintf('\nStep 7: Duffing alpha(T)\n');

%% ==== BEGIN PlotDuffingAlpha.m (integrated) ====

if exist('exportCache','var') && isfield(exportCache,'Monolithic') && ...
                                 isfield(exportCache,'Jointed')
    % Read from workspace cache — no file I/O needed
    obs_M = exportCache.Monolithic.obs;
    obs_J = exportCache.Jointed.obs;
    fprintf('PlotDuffingAlpha: reading from exportCache.\n');
else
    % Standalone: find the most recent exports automatically
    d_M = dir('export_Monolithic_observations_*.csv');
    d_J = dir('export_Jointed_observations_*.csv');
    if isempty(d_M) || isempty(d_J)
        error('PlotDuffingAlpha: no observation CSVs found. Run process2.m first.');
    end
    [~,iM] = max([d_M.datenum]);  [~,iJ] = max([d_J.datenum]);
    fprintf('PlotDuffingAlpha: reading %s\n', d_M(iM).name);
    fprintf('PlotDuffingAlpha: reading %s\n', d_J(iJ).name);
    obs_M = readtable(d_M(iM).name);
    obs_J = readtable(d_J(iJ).name);
end

% ── Compute alpha(T) for one table ────────────────────────────────────────

[T_M, alpha_M, ~, R2_M] = fitDuffingAlpha(obs_M);
[T_J, alpha_J, ~, R2_J] = fitDuffingAlpha(obs_J);

% ── Print summary ─────────────────────────────────────────────────────────
fprintf('\n=== Duffing alpha(T) — Monolithic ===\n');
fprintf('  %5s  %10s  %6s\n','T(C)','alpha','R2');
for i = 1:numel(T_M)
    if ~isnan(alpha_M(i))
        fprintf('  %5.0f  %+10.4f  %6.4f\n', T_M(i), alpha_M(i), R2_M(i));
    end
end

fprintf('\n=== Duffing alpha(T) — Jointed ===\n');
fprintf('  %5s  %10s  %6s\n','T(C)','alpha','R2');
for i = 1:numel(T_J)
    if ~isnan(alpha_J(i))
        fprintf('  %5.0f  %+10.4f  %6.4f\n', T_J(i), alpha_J(i), R2_J(i));
    end
end

% Interpolate zero-crossing of alpha_J
valid_J = ~isnan(alpha_J);
T_J_v   = T_J(valid_J);
a_J_v   = alpha_J(valid_J);
T_zero  = NaN;
for i = 1:numel(T_J_v)-1
    if a_J_v(i) < 0 && a_J_v(i+1) >= 0
        T_zero = T_J_v(i) + (-a_J_v(i)/(a_J_v(i+1)-a_J_v(i))) * ...
                 (T_J_v(i+1) - T_J_v(i));
        break;
    end
end
if ~isnan(T_zero)
    fprintf('\nJointed alpha zero-crossing: T ≈ %.1f°C\n', T_zero);
end

% ── Figure ────────────────────────────────────────────────────────────────
fig = figure('Color','w','Name','Duffing Coefficient vs Temperature', ...
    'Position',[100 100 700 450]);
ax = axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');

COL_M = [0.85 0.33 0.10];   % orange — monolithic
COL_J = [0.15 0.45 0.75];   % blue   — jointed

% Plot monolithic
valid_M = ~isnan(alpha_M);
plot(ax, T_M(valid_M), alpha_M(valid_M), '-o', ...
    'Color', COL_M, 'LineWidth', 2, 'MarkerSize', 7, ...
    'MarkerFaceColor', COL_M, 'DisplayName', 'Monolithic');

% Plot jointed — use open markers where R2 < 0.5 (two-mechanism regime, poor fit)
valid_J_good = ~isnan(alpha_J) & R2_J >= 0.5;
valid_J_poor = ~isnan(alpha_J) & R2_J <  0.5;

plot(ax, T_J(valid_J_good), alpha_J(valid_J_good), '-o', ...
    'Color', COL_J, 'LineWidth', 2, 'MarkerSize', 7, ...
    'MarkerFaceColor', COL_J, 'DisplayName', 'Jointed');

if any(valid_J_poor)
    plot(ax, T_J(valid_J_poor), alpha_J(valid_J_poor), 'o', ...
        'Color', COL_J, 'LineWidth', 2, 'MarkerSize', 7, ...
        'MarkerFaceColor', 'w', 'HandleVisibility','off');
    % Connect to the main line without legend entry
    all_J = ~isnan(alpha_J);
    plot(ax, T_J(all_J), alpha_J(all_J), '-', ...
        'Color', [COL_J, 0.4], 'LineWidth', 1, 'HandleVisibility','off');
end

% Reference lines
yline(ax, 0, 'k--', 'LineWidth', 1.2, 'HandleVisibility','off');
xline(ax, 0, 'k:', 'LineWidth', 1.2, 'HandleVisibility','off');

% Annotate zero-crossing
if ~isnan(T_zero)
    plot(ax, T_zero, 0, 'kd', 'MarkerSize', 10, 'MarkerFaceColor','k', ...
        'DisplayName', sprintf('\\alpha_J = 0  (T \\approx %.0f°C)', T_zero));
end

% Shade softening and hardening regions
xl = xlim(ax);
patch(ax, [xl(1) 0 0 xl(1)], [0 0 max(ylim(ax)) max(ylim(ax))], ...
    [0.9 0.75 0.75], 'FaceAlpha', 0.15, 'EdgeColor','none', ...
    'HandleVisibility','off');
patch(ax, [0 xl(2) xl(2) 0], [0 0 max(ylim(ax)) max(ylim(ax))], ...
    [0.75 0.85 0.95], 'FaceAlpha', 0.15, 'EdgeColor','none', ...
    'HandleVisibility','off');

% Region labels
text(ax, -25, max(ylim(ax))*0.85, 'Softening', ...
    'FontSize', 9, 'Color', [0.6 0.2 0.2], 'HorizontalAlignment','left');
text(ax,  5, max(ylim(ax))*0.85, 'Hardening', ...
    'FontSize', 9, 'Color', [0.2 0.4 0.7], 'HorizontalAlignment','left');

% Note: open markers = poor Duffing fit (two-mechanism regime)
if any(valid_J_poor)
    text(ax, min(T_J(valid_J_poor))-2, min(alpha_J(valid_J_poor))-0.3, ...
        'Open markers: R^2 < 0.5', 'FontSize', 8, 'Color', COL_J, ...
        'FontAngle','italic', 'HorizontalAlignment','center');
end

xlabel(ax, 'Temperature (°C)', 'FontSize', 12);
ylabel(ax, '\alpha  (Hz^2/\mum^2)', 'FontSize', 12);
title(ax, 'Temperature-Dependent Duffing Coefficient', 'FontSize', 13, ...
    'FontWeight', 'bold');
legend(ax, 'Location','northwest', 'FontSize', 10);
xlim(ax, [-35 125]);

% Save
outDir = fullfile('exports','figures');
if ~exist(outDir,'dir'); mkdir(outDir); end
saveas(fig, fullfile(outDir,'DuffingAlpha_vs_Temperature.png'));
saveas(fig, fullfile(outDir,'DuffingAlpha_vs_Temperature.fig'));
fprintf('Figure saved to %s\n', outDir);
%% ==== END PlotDuffingAlpha.m ====


%% ==== BEGIN couplingmodel.m (integrated) ====

T_test = [-30 -20 -10 0 10 20 30 40 50 60 70 80 90 100 110 120]';  % °C

alpha_M_meas = [-3.20 -2.70 -3.10 -3.40 -4.10 -4.77 -4.20 -3.80 ...
                -3.50 -3.20 -2.90 -2.60 -2.30 -2.00 -2.10 -2.20]';

alpha_J_meas = [-0.60  0.16  0.50  0.80  1.40  2.11  2.30  2.40 ...
                 2.49  2.35  2.10  1.90  1.80  1.70  1.60  2.20]';

R2_J = [0.12  0.03  0.73  0.78  0.82  0.88  0.90  0.91 ...
         0.92  0.90  0.88  0.85  0.83  0.82  0.80  0.79]';

%% ── 1. Beam and material parameters ──────────────────────────────────────

% Geometry (SI: metres, kg, N)
b   = 30e-3;        % beam width (m)
h   = 3e-3;         % beam thickness (m)
L_M = 330e-3;       % monolithic free span (m)
L_J = 360e-3;       % jointed effective span (m) — 2 x 180 mm sub-spans

A_beam = b * h;                     % cross-sectional area (m^2)
I_beam = b * h^3 / 12;              % second moment of area (m^4)

% Material (mild steel)
E       = 200e9;                    % Young's modulus (Pa)
rho     = 7850;                     % density (kg/m^3)
alpha_T = 12e-6;                    % CTE (1/°C)

% Bolt (M4 steel)
d_bolt  = 4e-3;                     % bolt shank diameter (m)
A_bolt  = pi * d_bolt^2 / 4;       % bolt shank area (m^2)
l_bolt  = 6e-3;                     % effective bolt grip length (m) — 2x3mm plates
n_bolts = 2;                        % number of bolts at joint

% Bolt and clamped-stack axial stiffness
k_bolt  = E * A_bolt / l_bolt;     % N/m
A_clamp = 3 * pi * (d_bolt/2)^2;  % effective clamped area (VDI 2230 cone approximation)
k_clamp = E * A_clamp / l_bolt;   % N/m

fprintf('Beam A = %.2e m^2,  I = %.2e m^4\n', A_beam, I_beam);
fprintf('k_bolt = %.1f kN/mm,  k_clamp = %.1f kN/mm\n', ...
        k_bolt/1e6, k_clamp/1e6);

%% ── 2. Identify thermal neutral point from monolithic data ────────────────
% T0 is where Delta_K(0.5V) changes sign.
% From paper: +0.238% at 20°C, -0.028% at 30°C → interpolate.

DeltaK_low_M = [0.238, -0.028];    % Delta_K(0.5V) at [20, 30] °C
T_bracket    = [20, 30];
T0 = interp1(DeltaK_low_M, T_bracket, 0);
fprintf('\nThermal neutral point T0 = %.1f °C\n', T0);

%% ── 3. Fit geometric softening model to monolithic data ───────────────────
% Model: alpha_geo(T) = alpha0 + beta * N(T)
%        N(T) = E * A * alpha_T * (T - T0)   [Pa * m^2 = N, but we work
%        in consistent Hz^2/um^2 units so N is normalised below]
%
% Rearranged for linear regression:
%        alpha_M = alpha0 + beta * [E*A*alpha_T*(T-T0)]
%
% N_vec: thermal axial force (N), positive = tensile

N_vec = E * A_beam * alpha_T * (T_test - T0);   % N

% Linear regression: [alpha0, beta]
X_geo  = [ones(size(T_test)), N_vec];
params_geo = X_geo \ alpha_M_meas;
alpha0_fit = params_geo(1);
beta_fit   = params_geo(2);

alpha_geo_fit = alpha0_fit + beta_fit * N_vec;

SS_res = sum((alpha_M_meas - alpha_geo_fit).^2);
SS_tot = sum((alpha_M_meas - mean(alpha_M_meas)).^2);
R2_geo = 1 - SS_res/SS_tot;

fprintf('\nGeometric softening fit:\n');
fprintf('  alpha0 = %.4f Hz^2/um^2\n', alpha0_fit);
fprintf('  beta   = %.4e Hz^2/um^2 per N\n', beta_fit);
fprintf('  R^2    = %.4f\n', R2_geo);

% Predicted transition from geo alone (where alpha_geo = 0)
if alpha0_fit * (alpha0_fit + beta_fit * E*A_beam*alpha_T*(200-T0)) < 0
    T_geo_cross = T0 - alpha0_fit / (beta_fit * E*A_beam*alpha_T);
    fprintf('  Geo-only zero-crossing: T = %.1f °C\n', T_geo_cross);
else
    T_geo_cross = NaN;
    fprintf('  Geo-only zero-crossing: not found in range\n');
end

%% ── 4. Extract contact hardening contribution ─────────────────────────────
% alpha_contact(T) = alpha_J(T) - alpha_geo(T)
% Only use temperatures where Duffing fit was reliable (R^2 > 0.7)

good_idx = R2_J > 0.70;
alpha_contact_meas = alpha_J_meas - alpha_geo_fit;

fprintf('\nContact hardening alpha_contact at selected temperatures:\n');
fprintf('  T(°C)   alpha_J   alpha_geo   alpha_contact   R^2_J\n');
for k = 1:numel(T_test)
    if good_idx(k)
        fprintf('  %4.0f    %+.3f     %+.3f       %+.3f           %.2f\n', ...
            T_test(k), alpha_J_meas(k), alpha_geo_fit(k), ...
            alpha_contact_meas(k), R2_J(k));
    end
end

%% ── 5. Model contact normal force P(T) ───────────────────────────────────
% P(T) = P0 + dP/dT * (T - T_assembly)
% For matched-material bolt and plates:
%   dP/dT = -(k_bolt * k_clamp)/(k_bolt + k_clamp) * A_bolt * E_bolt * alpha_T
%   (this is small for matched materials; P0 and T_assembly are free parameters)

dPdT = -(k_bolt * k_clamp)/(k_bolt + k_clamp) * alpha_T * E;   % N/°C per bolt
fprintf('\ndP/dT per bolt = %.2f N/°C\n', dPdT);

% P(T) as function of free parameters P0 (preload per bolt) and T_assembly
% P_fun: @(P0, T_a, T) P0 + n_bolts * dPdT * (T - T_a)
P_fun = @(P0, T_a, T) max(0, P0 + n_bolts * dPdT * (T - T_a));
% max(0,...) enforces that contact force cannot go negative (joint opens)

%% ── 6. Fit contact hardening model ───────────────────────────────────────
% Model: alpha_contact(T) = C * P(T)^gamma
% Free parameters: C, gamma, P0, T_assembly
% Fit only at temperatures with reliable Duffing fit (R^2_J > 0.7)

T_fit    = T_test(good_idx);
ac_fit   = alpha_contact_meas(good_idx);
w_fit    = R2_J(good_idx);              % weight by fit quality

% Only fit where alpha_contact > 0 (hardening regime)
pos_idx  = ac_fit > 0;
T_fit    = T_fit(pos_idx);
ac_fit   = ac_fit(pos_idx);
w_fit    = w_fit(pos_idx);

if numel(T_fit) < 3
    warning('CouplingModel: fewer than 3 reliable positive alpha_contact values — contact fit may be poorly conditioned.');
end

% Initial parameter guess
% P0: typical M4 bolt torqued to ~1 Nm → preload ~500 N
% T_assembly: room temperature 20°C
% gamma: ~1 (linear to start)
% C: alpha_contact / P0^gamma at reference temperature
p0_guess     = 500;
T_asm_guess  = 20;
gamma_guess  = 1.0;
C_guess      = mean(ac_fit) / p0_guess^gamma_guess;

x0 = [C_guess, gamma_guess, p0_guess, T_asm_guess];
lb = [0,    0.1,   10,   -50];
ub = [1e3,  5.0, 5000,    50];

obj_fun = @(x) sum(w_fit .* (ac_fit - x(1)*P_fun(x(3),x(4),T_fit).^x(2)).^2);

opts = optimset('Display','off','TolFun',1e-10,'TolX',1e-10,'MaxIter',5000);
try
    x_opt = fmincon(obj_fun, x0, [], [], [], [], lb, ub, [], ...
                    optimoptions('fmincon','Display','off','MaxIterations',5000, ...
                                 'OptimalityTolerance',1e-10,'StepTolerance',1e-10));
catch
    % Fallback to fminsearch if optimisation toolbox unavailable
    x_opt = fminsearch(obj_fun, x0, opts);
    x_opt(3) = max(10, x_opt(3));   % enforce P0 > 0
end

C_opt     = x_opt(1);
gamma_opt = x_opt(2);
P0_opt    = x_opt(3);
T_asm_opt = x_opt(4);

fprintf('\nContact hardening fit:\n');
fprintf('  C      = %.4e\n',  C_opt);
fprintf('  gamma  = %.4f\n',  gamma_opt);
fprintf('  P0     = %.1f N (per bolt: %.1f N)\n', P0_opt, P0_opt/n_bolts);
fprintf('  T_asm  = %.1f °C\n', T_asm_opt);

% Physical plausibility check for M4 bolt
P0_per_bolt = P0_opt / n_bolts;
if P0_per_bolt > 100 && P0_per_bolt < 3000
    fprintf('  P0 per bolt within plausible range for M4 bolt.\n');
else
    warning('CouplingModel: fitted P0 (%.0f N/bolt) outside expected range for M4 — check input data.', P0_per_bolt);
end

%% ── 7. Construct full model over continuous temperature range ─────────────

T_cont = linspace(-40, 130, 500)';

% Geometric softening (from monolithic fit)
N_cont         = E * A_beam * alpha_T * (T_cont - T0);
alpha_geo_cont = alpha0_fit + beta_fit * N_cont;

% Contact hardening (from fitted model)
P_cont             = P_fun(P0_opt, T_asm_opt, T_cont);
alpha_contact_cont = C_opt * P_cont.^gamma_opt;

% Full superposition
alpha_J_model = alpha_geo_cont + alpha_contact_cont;

%% ── 8. Find predicted transition temperatures ─────────────────────────────

% Geo-only transition
if ~isnan(T_geo_cross)
    fprintf('\nPredicted transition temperatures:\n');
    fprintf('  Geometric softening alone:    T = %.1f °C\n', T_geo_cross);
else
    fprintf('\n  Geometric softening alone:    no zero-crossing in range\n');
end

% Full model transition
sign_changes = find(diff(sign(alpha_J_model)));
if ~isempty(sign_changes)
    T_model_cross = interp1(alpha_J_model(sign_changes(1):sign_changes(1)+1), ...
                            T_cont(sign_changes(1):sign_changes(1)+1), 0);
    fprintf('  Full superposition model:     T = %.1f °C\n', T_model_cross);
    fprintf('  Observed (BFS sign change):   T = 0 °C\n');
    fprintf('  Observed (Duffing crossing):  T ≈ -22 °C\n');
    fprintf('  Gap geo-only vs observed BFS: %.1f °C\n', T_geo_cross - 0);
    fprintf('  Gap model  vs observed BFS:   %.1f °C\n', T_model_cross - 0);
else
    fprintf('  Full model: no zero-crossing found — check input Duffing data.\n');
end

%% ── 9. Plot ───────────────────────────────────────────────────────────────

figure('Name','Coupling Model — Duffing Coefficient vs Temperature', ...
       'Position',[100 100 820 520]);

% Shaded background: softening (left) vs hardening (right)
hold on;
fill([T_cont(1) T_cont(end) T_cont(end) T_cont(1)], [0 0 6 6], ...
     [0.9 1.0 0.9], 'EdgeColor','none', 'FaceAlpha',0.4);
fill([T_cont(1) T_cont(end) T_cont(end) T_cont(1)], [-6 -6 0 0], ...
     [1.0 0.9 0.9], 'EdgeColor','none', 'FaceAlpha',0.4);
yline(0, 'k-', 'LineWidth', 0.8);

% Three model curves
p1 = plot(T_cont, alpha_geo_cont,   'b--', 'LineWidth', 1.8);
p2 = plot(T_cont, alpha_contact_cont,'r:', 'LineWidth', 1.8);
p3 = plot(T_cont, alpha_J_model,    'k-', 'LineWidth', 2.2);

% Measured data points
p4 = scatter(T_test, alpha_M_meas, 40, 'b', 'filled', 'DisplayName', ...
             'Monolithic measured $\alpha_M$');
p5 = scatter(T_test(good_idx), alpha_J_meas(good_idx), 40, 'k', 'filled', ...
             'DisplayName', 'Jointed measured $\alpha_J$ (reliable fit)');
scatter(T_test(~good_idx), alpha_J_meas(~good_idx), 40, 'k', ...
        'LineWidth', 1.2, 'DisplayName', 'Jointed measured $\alpha_J$ (poor fit)');

% Transition temperature markers
if ~isnan(T_geo_cross)
    xline(T_geo_cross, 'b--', sprintf('Geo only\n%.0f°C', T_geo_cross), ...
          'LineWidth', 1.2, 'LabelVerticalAlignment','bottom', ...
          'LabelHorizontalAlignment','left');
end
if ~isempty(sign_changes)
    xline(T_model_cross, 'k-', sprintf('Model\n%.0f°C', T_model_cross), ...
          'LineWidth', 1.2, 'LabelVerticalAlignment','bottom');
end
xline(0,   'k:', 'BFS sign change (0°C)',   'LineWidth', 1.0, ...
      'LabelVerticalAlignment','top');
xline(-22, 'k:', '\alpha_J cross (–22°C)', 'LineWidth', 1.0, ...
      'LabelVerticalAlignment','top', 'LabelHorizontalAlignment','right');

% Labels and formatting
xlabel('Temperature (°C)', 'FontSize', 12);
ylabel('Duffing coefficient $\alpha$ (Hz$^2$/$\mu$m$^2$)', ...
       'Interpreter','latex', 'FontSize', 12);
title('Decomposition of Duffing coefficient: geometric softening vs contact hardening', ...
      'FontSize', 11);
% Use MATLAB's supported TeX markup; avoid LaTeX math delimiters and
% Unicode punctuation in the figure legend (vector rendering safe).
legend([p1 p2 p3 p4 p5], ...
       {'\alpha_{geo}(T) - geometric softening (monolithic fit)', ...
        '\alpha_{contact}(T) - contact hardening (Mindlin fit)', ...
        '\alpha_J(T) = \alpha_{geo} + \alpha_{contact} - full model', ...
        'Monolithic measured', ...
        'Jointed measured (reliable)'}, ...
       'Interpreter','tex', 'Location','southwest', 'FontSize', 9);
text(-38, 4.5, 'Hardening', 'Color',[0.2 0.6 0.2], 'FontSize',10);
text(-38,-4.5, 'Softening', 'Color',[0.7 0.2 0.2], 'FontSize',10);
xlim([-40 130]);
grid on; box on;

% Use an absolute path for MATLAB's PDF output device on Windows.
% If vector export is unavailable, save a high-resolution PNG instead.
figureOutputDir = fullfile(pwd, 'exports', 'figures');
if ~isfolder(figureOutputDir)
    [directoryCreated, directoryMessage] = mkdir(figureOutputDir);
    if ~directoryCreated
        error('research_pipeline:OutputDirectory', ...
              'Cannot create output directory %s: %s', figureOutputDir, directoryMessage);
    end
end
figurePdfPath = fullfile(figureOutputDir, 'CouplingModel_DuffingDecomposition.pdf');
try
    exportgraphics(gcf, figurePdfPath, 'ContentType', 'vector');
    fprintf('Figure exported: %s\n', figurePdfPath);
catch pdfExportError
    warning('research_pipeline:PDFExport', ...
            'PDF export failed (%s). Attempting PNG fallback.', pdfExportError.message);
    figurePngPath = fullfile(figureOutputDir, 'CouplingModel_DuffingDecomposition.png');
    exportgraphics(gcf, figurePngPath, 'Resolution', 300);
    fprintf('Figure exported as PNG: %s\n', figurePngPath);
end

%% ── 10. Summary table ─────────────────────────────────────────────────────

fprintf('\n══════════════════════════════════════════════════════════\n');
fprintf('  Coupling model summary\n');
fprintf('══════════════════════════════════════════════════════════\n');
fprintf('  Monolithic neutral point (from Delta_K data):  %.1f °C\n', T0);
fprintf('  Geo-only predicted transition:                 %.1f °C\n', T_geo_cross);
if ~isempty(sign_changes)
fprintf('  Two-mechanism model predicted transition:      %.1f °C\n', T_model_cross);
end
fprintf('  Observed BFS sign reversal:                    0.0 °C\n');
fprintf('  Observed Duffing zero-crossing:               -22.0 °C\n');
fprintf('\n  Gap (geo-only vs BFS):     %.1f °C\n', T_geo_cross - 0);
if ~isempty(sign_changes)
fprintf('  Gap (model vs BFS):        %.1f °C  ', T_model_cross - 0);
if abs(T_model_cross) < 10
    fprintf('[within experimental range]\n');
else
    fprintf('[outside experimental range — review input data]\n');
end
end
fprintf('══════════════════════════════════════════════════════════\n');
%% ==== END couplingmodel.m ====

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ── DeltaK figure — both systems overlaid ─────────────────────────────────
rep_temps = [-30 20 70 120];
T_vec_plot = config.steppedsine.temperature.Vector;

figure('Name', 'Amplitude-Dependent Stiffness Deviation');
for sp = 1:2
    subplot(1,2,sp); hold on; grid on; box on;
    if sp == 1
        sysLabel = 'Monolithic'; mkr = 'o';
    else
        sysLabel = 'Jointed';    mkr = 's';
    end
    xlabel('Displacement Amplitude (\mum)');
    ylabel('\DeltaK / K_0 (%)');
    title(sysLabel);
    yline(0, '--', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');

    for t = 1:length(rep_temps)
        t_idx = find(T_vec_plot == rep_temps(t));
        if isempty(t_idx); continue; end

        A_um  = NaN(1, config.steppedsine.nTemps);   % placeholder size
        dK_pct = NaN(1, 16);
        bb = proc_all.(sysLabel).steppedsine;

        % Rebuild per-amplitude deltaK and displacement from stored amp structs
        Na_loc = numel(bb.amp{t_idx});
        A_vec  = NaN(1, Na_loc);
        dK_vec = NaN(1, Na_loc);
        for a = 1:Na_loc
            amp_a   = bb.amp{t_idx}(a);
            dK_vec(a) = amp_a.Stiffness.deltaK * 100;   % convert to %
            % Displacement from backbone A_pts (metres -> µm)
            A_vec(a)  = proc_all.(sysLabel).steppedsine.backbone{t_idx}.A_pts(a) * 1e6;
        end

        valid = ~isnan(A_vec) & ~isnan(dK_vec) & A_vec > 0;
        if sum(valid) < 2; continue; end

        plot(A_vec(valid), dK_vec(valid), [mkr '-'], ...
            'LineWidth', 1.5, 'MarkerSize', 4, ...
            'DisplayName', sprintf('%d°C', rep_temps(t)));
    end
    legend('Location', 'best');
end

exportgraphics(gcf, fullfile('exports','figures','deltaK_combined.png'), 'Resolution', 300);

% ── BFS vs Temperature — Combined Overlay ─────────────────────────────────
fprintf('\nBFS overlay plot\n');

T_vec_plot = config.steppedsine.temperature.Vector;

% Extract BFS for each system from stored proc_all
systems = {'Monolithic', 'Jointed'};
colors  = {[0 0 0], [0.13 0.47 0.71]};   % black, blue — matches Figure 2
markers = {'o', 's'};

figure('Name', 'BFS vs Temperature — Combined');
hold on; grid on; box on;

yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, 'HandleVisibility', 'off');

for s = 1:2
    sysLabel = systems{s};
    bb       = proc_all.(sysLabel).steppedsine.backbone;
    BFS_vec  = NaN(1, config.steppedsine.nTemps);

    for Z = 1:config.steppedsine.nTemps
        f_pts   = bb{Z}.f_pts;
        valid_f = ~isnan(f_pts);
        if sum(valid_f) < 2; continue; end
        f_valid    = f_pts(valid_f);
        BFS_vec(Z) = (f_valid(end) - f_valid(1)) / f_valid(1) * 100;
    end

    proc_all.(sysLabel).thermo.BFS.vec = BFS_vec;
    proc_all.(sysLabel).thermo.BFS.T   = T_vec_plot;

    plot(T_vec_plot, BFS_vec, ...
        [markers{s} '-'],            ...
        'Color',           colors{s}, ...
        'LineWidth',       1.5,        ...
        'MarkerFaceColor', colors{s}, ...
        'MarkerSize',      5,          ...
        'DisplayName',     sysLabel);
end

% Annotate sign change on jointed curve only
BFS_joint = proc_all.Jointed.thermo.BFS.vec;
sc_idx    = find(diff(sign(BFS_joint)) ~= 0, 1);
if ~isempty(sc_idx)
    xc = mean(T_vec_plot(sc_idx:sc_idx+1));
    xline(xc, ':', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.2, ...
        'Label',                    sprintf('Sign change \\approx %.0f°C', xc), ...
        'LabelVerticalAlignment',   'bottom', ...
        'HandleVisibility',         'off');
end

xlabel('Temperature (°C)', 'FontSize', 11);
ylabel('BFS (%)',           'FontSize', 11);
xlim([min(T_vec_plot)-5, max(T_vec_plot)+5]);
legend('Location', 'northwest');

if ~exist(fullfile('exports','figures'), 'dir')
    mkdir(fullfile('exports','figures'));
end
exportgraphics(gcf, fullfile('exports','figures','BFS_vs_T_combined.png'), ...
    'Resolution', 300);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all

% run PaperFigures.m
COL_M = [0    0    0   ];   % monolithic = black
COL_J = [0    0    1   ];   % jointed    = blue
% =========================================================================
% OVERLAY BACKBONE PLOTS — 4 temperatures on one figure, per system
% =========================================================================

rep_temps  = [-30, 20, 70, 120];
rep_colors = lines(numel(rep_temps));
T_vec      = config.steppedsine.temperature.Vector;

for sp = 1:2
    if sp == 1; sysLabel = 'Monolithic'; else; sysLabel = 'Jointed'; end

    % Recover the freq range from the filenames for xlim
    if strcmp(sysLabel, 'Monolithic')
        fStart = 147; fEnd = 152; %proc_all.Monolithic.impact.monolithic.peakFreq;
    else
        fStart = 112; fEnd = 117; %proc_all.Jointed.impact.jointed.peakFreq;
    end

    figOV = figure('Name', sprintf('Backbone Overlay — %s', sysLabel));
    ax_ov = axes(figOV);
    hold(ax_ov, 'on');
    grid(ax_ov, 'on');
    set(ax_ov, 'YScale', 'log');
    % xlim(ax_ov, [0.997 1.003]);   % adjust range to suit your data
    ylim(ax_ov, [0.1 100]);
    xlabel(ax_ov, 'Normalised frequency  f / f_0', 'FontSize', 11);
    ylabel(ax_ov, 'Response Amplitude (g)', 'FontSize', 11);

    for ti = 1:numel(rep_temps)
    T_r = rep_temps(ti);
    Z_r = find(T_vec == T_r, 1);
    if isempty(Z_r); continue; end

    col  = rep_colors(ti,:);
    Na_r = numel(proc_all.(sysLabel).steppedsine.amp{Z_r});
    bb_r = proc_all.(sysLabel).steppedsine.backbone{Z_r};

    % f0 = low-amplitude resonance frequency at this temperature
    % Use the first valid backbone point as the reference
    f_pts_r = bb_r.f_pts;
    valid_f = ~isnan(f_pts_r);
    if sum(valid_f) < 1; continue; end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f0_sq = proc_all.(sysLabel).thermo.duffing.f0sq(Z_r);
if isnan(f0_sq) || f0_sq <= 0
    f0 = fStart+(fEnd-fStart)/2;  % fallback
else
    f0 = sqrt(f0_sq);
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % FRF curves — normalise x by f0
    for a = 1:Na_r
        amp_a = proc_all.(sysLabel).steppedsine.amp{Z_r}(a);
        if ~isempty(amp_a.Fit.fDense) && ~isempty(amp_a.Fit.yFit)
            plot(ax_ov, amp_a.Fit.fDense / f0, amp_a.Fit.yFit, '-', ...
                'Color', [col 0.20], 'LineWidth', 0.8, ...
                'HandleVisibility', 'off');
        end
    end

    % Backbone curve — normalise x by f0
    if ~all(isnan(bb_r.f)) && ~all(isnan(bb_r.y))
        plot(ax_ov, bb_r.f / f0, bb_r.y, '-', ...
            'Color', col, 'LineWidth', 2, ...
            'DisplayName', sprintf('%d°C', T_r));

        % Duffing extrapolation to zero amplitude
        alpha_T = proc_all.(sysLabel).thermo.duffing.alpha(Z_r);
        if ~isnan(alpha_T)
            valid_bb = ~isnan(bb_r.A_pts) & ~isnan(bb_r.y_pts);
            if sum(valid_bb) >= 1
                A_min_um = bb_r.A_pts(find(valid_bb, 1, 'first')) * 1e6;
                y_min    = bb_r.y_pts(find(valid_bb, 1, 'first'));
                A_extrap = linspace(0, A_min_um, 30);
                f_extrap = sqrt(f0_sq + alpha_T .* A_extrap.^2);
                y_extrap = y_min .* (A_extrap / A_min_um);
                y_extrap(1) = NaN;
                plot(ax_ov, f_extrap / f0, y_extrap, '--', ...
                    'Color', col, 'LineWidth', 1.2, ...
                    'HandleVisibility', 'off');
            end
        end
    end
    end

    legend(ax_ov, 'Location', 'northwest', 'FontSize', 10);

    exportgraphics(figOV, ...
        fullfile('exports', 'figures', 'paper', ...
            sprintf('backbone_overlay_%s.png', sysLabel)), ...
        'Resolution', 300);
    fprintf('Backbone overlay saved: %s\n', sysLabel);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% =========================================================================
% BACKBONE OVERLAY — Jointed vs Monolithic at 20°C
% =========================================================================

T_target  = 20;
Z_20_M    = find(T_vec == T_target, 1);
Z_20_J    = find(T_vec == T_target, 1);

figComp = figure('Name', 'Backbone Overlay — Jointed vs Monolithic at 20°C');
ax_comp = axes(figComp);
hold(ax_comp, 'on');
grid(ax_comp, 'on');
set(ax_comp, 'YScale', 'log');
% ylim(ax_comp, [0 100]);
xlabel(ax_comp, 'Normalised frequency  f / f_0', 'FontSize', 11);
ylabel(ax_comp, 'Response Amplitude (g)', 'FontSize', 11);

systems_comp = {'Monolithic', 'Jointed'};
colors_comp  = {COL_M, COL_J};
% fStarts      = [147, 112];
% fEnds        = [152, 117];
Z_comp       = {Z_20_M, Z_20_J};

for sp = 1:2

    if sp == 1; sysLabel = 'Monolithic'; else; sysLabel = 'Jointed'; end

    col      = colors_comp{sp};
    Z_r      = Z_comp{sp};

    Na_r = numel(proc_all.(sysLabel).steppedsine.amp{Z_r});
    bb_r = proc_all.(sysLabel).steppedsine.backbone{Z_r};

    % Recover the freq range from the filenames for xlim
    if strcmp(sysLabel, 'Monolithic')
        fStart = 147; fEnd = 152; %proc_all.Monolithic.impact.monolithic.peakFreq;
    else
        fStart = 112; fEnd = 117; %proc_all.Jointed.impact.jointed.peakFreq;
    end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f0_sq = proc_all.(sysLabel).thermo.duffing.f0sq(Z_r);
if isnan(f0_sq) || f0_sq <= 0
    f0 = fStart+(fEnd-fStart)/2;  % fallback
else
    f0 = sqrt(f0_sq);
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % FRF curves — light version of system colour
    for a = 1:Na_r
        amp_a = proc_all.(sysLabel).steppedsine.amp{Z_r}(a);
        if ~isempty(amp_a.Fit.fDense) && ~isempty(amp_a.Fit.yFit)
            plot(ax_comp, amp_a.Fit.fDense  / f0, amp_a.Fit.yFit, '-', ...
                'Color', [col 0.20], 'LineWidth', 0.8, ...
                'HandleVisibility', 'off');
        end
    end

    % Backbone curve — solid, labelled
    if ~all(isnan(bb_r.f)) && ~all(isnan(bb_r.y))
        plot(ax_comp, bb_r.f / f0, bb_r.y, '-', ...
            'Color', col, 'LineWidth', 2, ...
            'DisplayName', sysLabel);

        % Extend backbone to f/f0 = 1 using lowest amplitude FRF peak
        alpha_T = proc_all.(sysLabel).thermo.duffing.alpha(Z_r);
        if ~isnan(alpha_T)
            % Find the lowest amplitude FRF peak response at this temperature
            y_lin = NaN;
            for a_ext = 1:Na_r
                amp_a = proc_all.(sysLabel).steppedsine.amp{Z_r}(a_ext);
                if ~isempty(amp_a.Fit.yFit) && ~all(isnan(amp_a.Fit.yFit))
                    y_lin = max(amp_a.Fit.yFit);  % peak of lowest amplitude FRF
                    break
                end
            end

            if ~isnan(y_lin) && ~isnan(f0_sq)
                % Get lowest backbone point
                valid_bb = ~isnan(bb_r.A_pts) & ~isnan(bb_r.y_pts);
                if sum(valid_bb) >= 1
                    A_min_um = bb_r.A_pts(find(valid_bb, 1, 'first')) * 1e6;
                    y_min    = bb_r.y_pts(find(valid_bb, 1, 'first'));
                    % Extrapolate frequency from Duffing fit
                    A_extrap = linspace(0, A_min_um, 30);
                    f_extrap = sqrt(f0_sq + alpha_T .* A_extrap.^2);
                    % Interpolate response between y_lin (at A=0) and y_min (at A_min)
                    y_extrap = linspace(y_lin, y_min, 30);
                    plot(ax_ov, f_extrap / f0, y_extrap, '--', ...
                        'Color', col, 'LineWidth', 1.2, ...
                        'HandleVisibility', 'off');
                end
            end
        end
    end
end

% xlim is tricky here since the two systems span different frequency ranges
% — leave it auto so both are visible, or set manually once you see the output
% xlim(ax_comp, 'auto');

legend(ax_comp, 'Location', 'northwest', 'FontSize', 10);

exportgraphics(figComp, ...
    fullfile('exports', 'figures', 'paper', ...
        'backbone_overlay_JointedVsMono_20C.png'), ...
    'Resolution', 300);
fprintf('Backbone overlay Jointed vs Monolithic at 20°C saved.\n');

%% Generate Q&A backup plots
close all

%% ==== BEGIN QA_BackupPlots.m (integrated) ====
%% QA_BackupPlots.m
% =========================================================================
% Q&A / BACKUP FIGURES FOR ISMA 2026
%
% Run AFTER process2.m has generated:
%   proc_all.Monolithic
%   proc_all.Jointed
%   config.steppedsine.temperature.Vector
%
% Produces:
%   1) Duffing R^2 vs temperature
%   2) BFS vs temperature + endpoint processing illustration
%   3) PCHIP peak-interpolation sensitivity example
%
% Output:
%   exports/figures/QA/
% =========================================================================

fprintf('\n============================================================\n');
fprintf('Generating Q&A backup figures...\n');
fprintf('============================================================\n');

qaDir = fullfile('exports','figures','QA');
if ~exist(qaDir,'dir')
    mkdir(qaDir);
end

systems = {'Monolithic','Jointed'};

% Prefer temperature vector stored with Duffing results.
if isfield(proc_all.Jointed,'thermo') && ...
   isfield(proc_all.Jointed.thermo,'duffing') && ...
   isfield(proc_all.Jointed.thermo.duffing,'T_vec')

    T = proc_all.Jointed.thermo.duffing.T_vec(:);

else
    T = config.steppedsine.temperature.Vector(:);
end

%% =========================================================================
% 1. DUFFING FIT QUALITY (R^2) VS TEMPERATURE
% =========================================================================

figR2 = figure('Name','Q&A - Duffing Fit Quality',...
               'Color','w');

ax = axes(figR2);
hold(ax,'on');
grid(ax,'on');
box(ax,'on');

for s = 1:numel(systems)

    sys = systems{s};

    if ~isfield(proc_all.(sys),'thermo') || ...
       ~isfield(proc_all.(sys).thermo,'duffing') || ...
       ~isfield(proc_all.(sys).thermo.duffing,'R2')

        warning('%s Duffing R2 data not found.',sys);
        continue;
    end

    R2 = proc_all.(sys).thermo.duffing.R2(:);

    if strcmpi(sys,'Monolithic')
        plot(ax,T,R2,'o-',...
            'LineWidth',1.7,...
            'MarkerSize',6,...
            'MarkerFaceColor','auto',...
            'DisplayName','Monolithic');
    else
        plot(ax,T,R2,'s-',...
            'LineWidth',1.7,...
            'MarkerSize',6,...
            'MarkerFaceColor','auto',...
            'DisplayName','Jointed');
    end
end

yline(ax,1,'--','HandleVisibility','off');

% Useful guide rather than a formal statistical threshold.
yline(ax,0.7,':','R^2 = 0.7 guide',...
    'LabelHorizontalAlignment','left',...
    'HandleVisibility','off');

xline(ax,0,':','0 ^\circC',...
    'LabelVerticalAlignment','bottom',...
    'HandleVisibility','off');

xlabel(ax,'Temperature (^oC)');
ylabel(ax,'Duffing fit R^2');

title(ax,'Cubic Duffing Fit Quality vs Temperature');

xlim(ax,[min(T)-5 max(T)+5]);
ylim(ax,[0 1.05]);

legend(ax,'Location','best');

% Print exact values to command window.
fprintf('\nDUFFING FIT QUALITY\n');

for s = 1:numel(systems)

    sys = systems{s};

    if isfield(proc_all.(sys).thermo.duffing,'R2')

        R2 = proc_all.(sys).thermo.duffing.R2(:);

        fprintf('\n%s\n',sys);
        fprintf('Temperature [degC]      R2\n');

        for i = 1:min(numel(T),numel(R2))
            fprintf('%8.0f              %.4f\n',T(i),R2(i));
        end
    end
end

exportgraphics(figR2,...
    fullfile(qaDir,'QA_Duffing_R2_vs_T.png'),...
    'Resolution',300);


%% =========================================================================
% 2. BFS — DIRECTLY FROM BACKBONE ENDPOINTS
% =========================================================================

figBFS = figure('Name','Q&A - BFS Definition',...
                'Color','w');

ax = axes(figBFS);
hold(ax,'on');
grid(ax,'on');
box(ax,'on');

BFS_store = struct();

for s = 1:numel(systems)

    sys = systems{s};

    bb = proc_all.(sys).steppedsine.backbone;

    BFS = NaN(size(T));

    for Z = 1:min(numel(T),numel(bb))

        f_pts = bb{Z}.f_pts(:);

        valid = ~isnan(f_pts);

        if sum(valid) < 2
            continue
        end

        f_valid = f_pts(valid);

        % -------------------------------------------------------------
        % EXACT BFS DEFINITION USED IN process2.m
        % -------------------------------------------------------------
        f_lowAmp  = f_valid(1);
        f_highAmp = f_valid(end);

        BFS(Z) = ...
            (f_highAmp - f_lowAmp) / f_lowAmp * 100;

    end

    BFS_store.(sys) = BFS;

    if strcmpi(sys,'Monolithic')
        plot(ax,T,BFS,'o-',...
            'LineWidth',1.7,...
            'MarkerSize',6,...
            'MarkerFaceColor','auto',...
            'DisplayName','Monolithic');
    else
        plot(ax,T,BFS,'s-',...
            'LineWidth',1.7,...
            'MarkerSize',6,...
            'MarkerFaceColor','auto',...
            'DisplayName','Jointed');
    end
end

yline(ax,0,'--',...
    'HandleVisibility','off');

xlabel(ax,'Temperature (^oC)');
ylabel(ax,'BFS (%)');

title(ax,...
    'Backbone Frequency Shift: Endpoint Metric');

legend(ax,'Location','best');

xlim(ax,[min(T)-5 max(T)+5]);

% Find jointed zero crossing interval.
if isfield(BFS_store,'Jointed')

    bj = BFS_store.Jointed;

    idx = find(...
        ~isnan(bj(1:end-1)) & ...
        ~isnan(bj(2:end)) & ...
        bj(1:end-1).*bj(2:end) <= 0,...
        1,'first');

    if ~isempty(idx)

        T1 = T(idx);
        T2 = T(idx+1);

        fprintf('\nJointed BFS changes sign between %.0f and %.0f degC.\n',...
            T1,T2);

        xline(ax,mean([T1 T2]),':',...
            sprintf('sign reversal: %.0f to %.0f ^oC',T1,T2),...
            'LabelVerticalAlignment','bottom',...
            'HandleVisibility','off');
    end
end

fprintf('\nBFS DEFINITION\n');
fprintf('BFS = (f_highAmp - f_lowAmp) / f_lowAmp * 100 %%\n');

exportgraphics(figBFS,...
    fullfile(qaDir,'QA_BFS_vs_T.png'),...
    'Resolution',300);


%% =========================================================================
% 2B. BFS PROCESSING ILLUSTRATION AT ONE REPRESENTATIVE TEMPERATURE
%
% Shows:
%   measured backbone points
%   first valid point
%   last valid point
%   the exact quantities entering BFS
%
% Change this if you want another example temperature.
% =========================================================================

T_demo = 20;

[~,Zdemo] = min(abs(T-T_demo));
T_demo_actual = T(Zdemo);

figBFSDemo = figure('Name','Q&A - BFS Endpoint Definition',...
                    'Color','w');

ax = axes(figBFSDemo);
hold(ax,'on');
grid(ax,'on');
box(ax,'on');

sys = 'Jointed';

bb = proc_all.(sys).steppedsine.backbone{Zdemo};

fpts = bb.f_pts(:);

if isfield(bb,'A_pts')
    Aum = bb.A_pts(:)*1e6;
else
    Aum = (1:numel(fpts))';
end

valid = ~isnan(fpts) & ~isnan(Aum);

fvalid = fpts(valid);
Avalid = Aum(valid);

plot(ax,Avalid,fvalid,'o-',...
    'LineWidth',1.5,...
    'MarkerSize',6,...
    'DisplayName','Measured resonance points');

if numel(fvalid) >= 2

    fLow  = fvalid(1);
    fHigh = fvalid(end);

    ALow  = Avalid(1);
    AHigh = Avalid(end);

    BFS_demo = (fHigh-fLow)/fLow*100;

    scatter(ax,ALow,fLow,100,'filled',...
        'DisplayName','Lowest-amplitude point');

    scatter(ax,AHigh,fHigh,100,'filled',...
        'DisplayName','Highest-amplitude point');

    text(ax,ALow,fLow,...
        sprintf('  f_{low} = %.3f Hz',fLow),...
        'VerticalAlignment','bottom');

    text(ax,AHigh,fHigh,...
        sprintf('  f_{high} = %.3f Hz',fHigh),...
        'VerticalAlignment','top');

    title(ax,sprintf(...
        'BFS Endpoint Calculation — Jointed, %.0f ^oC | BFS = %.3f%%',...
        T_demo_actual,BFS_demo));
else

    title(ax,sprintf(...
        'BFS Endpoint Calculation — insufficient points at %.0f ^oC',...
        T_demo_actual));
end

xlabel(ax,'Response displacement amplitude (\mum)');
ylabel(ax,'Resonance frequency (Hz)');

legend(ax,'Location','best');

exportgraphics(figBFSDemo,...
    fullfile(qaDir,...
    sprintf('QA_BFS_endpoint_demo_T%+03.0f.png',T_demo_actual)),...
    'Resolution',300);


%% =========================================================================
% 3. PCHIP INTERPOLATION SENSITIVITY
%
% This plot compares:
%   a) raw 0.1-Hz sampled response
%   b) PCHIP interpolation used in SteppedSine2.m
%   c) raw-grid maximum
%   d) interpolated maximum
%
% It answers:
% "Could interpolation artificially shift the calculated resonance peak?"
%
% The difference is printed explicitly.
% =========================================================================

T_interp_demo = 20;

[~,Zi] = min(abs(T-T_interp_demo));
T_interp_actual = T(Zi);

sys = 'Jointed';

ampStruct = proc_all.(sys).steppedsine.amp{Zi};

% Use highest valid amplitude by default because nonlinear peak curvature
% is usually easiest to see there.
aDemo = numel(ampStruct);

while aDemo >= 1 && ...
     (isempty(ampStruct(aDemo).Fit.fDense) || ...
      isempty(ampStruct(aDemo).Fit.yFit))

    aDemo = aDemo - 1;
end

if aDemo >= 1

    % Stored PCHIP curve from SteppedSine2.m
    fDense = ampStruct(aDemo).Fit.fDense(:);
    yDense = ampStruct(aDemo).Fit.yFit(:);

    [yPchipMax,iPchip] = max(yDense);
    fPchip = fDense(iPchip);

    % -------------------------------------------------------------
    % Reconstruct the original ~0.1 Hz samples.
    %
    % proc_all stores the smoothed / processed FRF structures, but the
    % resonance-picking routine works from proc.peaks.ResponseAcc.
    %
    % If raw ResponseAcc was NOT copied into proc_all, infer the original
    % frequency nodes by sampling the stored PCHIP curve at the configured
    % stepped-sine frequency grid.
    %
    % This does NOT independently validate PCHIP; it is a visual sensitivity
    % diagnostic. See note below.
    % -------------------------------------------------------------

    freqStep = config.steppedsine.freq.Step;

    fRaw = ampStruct(aDemo).Fit.fRaw;
yRaw = ampStruct(aDemo).Fit.yRaw;

valid = ~isnan(fRaw) & ~isnan(yRaw);
fRaw = fRaw(valid);
yRaw = yRaw(valid);

% Maximum available directly from the measured 0.1 Hz grid
[yRawMax, iRaw] = max(yRaw);
fRawPeak = fRaw(iRaw);

% PCHIP result already used by your processing
fDense = ampStruct(aDemo).Fit.fDense;
yDense = ampStruct(aDemo).Fit.yFit;

[yInterpMax, iInterp] = max(yDense);
fInterpPeak = fDense(iInterp);

fprintf('Raw-grid peak = %.4f Hz\n', fRawPeak);
fprintf('PCHIP peak    = %.4f Hz\n', fInterpPeak);
fprintf('Shift         = %.4f Hz\n', fInterpPeak-fRawPeak);

    [yRawMax,iRaw] = max(yRaw);
    fRawPeak = fRaw(iRaw);

    df = fPchip - fRawPeak;

    figInterp = figure(...
        'Name','Q&A - PCHIP Peak Sensitivity',...
        'Color','w');

    ax = axes(figInterp);
    hold(ax,'on');
    grid(ax,'on');
    box(ax,'on');

    plot(ax,fRaw,yRaw,'o',...
        'MarkerSize',6,...
        'DisplayName','0.1 Hz frequency nodes');

    plot(ax,fDense,yDense,'-',...
        'LineWidth',1.8,...
        'DisplayName','PCHIP interpolation');

    scatter(ax,fRawPeak,yRawMax,100,'filled',...
        'DisplayName',sprintf(...
        'Discrete-grid max = %.3f Hz',fRawPeak));

    scatter(ax,fPchip,yPchipMax,100,'filled',...
        'DisplayName',sprintf(...
        'PCHIP max = %.3f Hz',fPchip));

    xlabel(ax,'Frequency (Hz)');
    ylabel(ax,'Response amplitude (g)');

    title(ax,sprintf(...
        ['Interpolation sensitivity — %s, %.0f ^oC, amplitude level %d\n' ...
         '\\Deltaf_{peak} = %.4f Hz'],...
        sys,T_interp_actual,aDemo,df));

    legend(ax,'Location','best');

    fprintf('\nPCHIP PEAK SENSITIVITY EXAMPLE\n');
    fprintf('System               : %s\n',sys);
    fprintf('Temperature          : %.0f degC\n',T_interp_actual);
    fprintf('Amplitude index      : %d\n',aDemo);
    fprintf('Discrete-grid peak   : %.6f Hz\n',fRawPeak);
    fprintf('PCHIP peak           : %.6f Hz\n',fPchip);
    fprintf('Difference            : %+.6f Hz\n',df);
    fprintf('Raw frequency spacing : %.6f Hz\n',freqStep);

    exportgraphics(figInterp,...
        fullfile(qaDir,...
        sprintf('QA_PCHIP_peak_sensitivity_%s_T%+03.0f.png',...
        sys,T_interp_actual)),...
        'Resolution',300);

else
    warning('Could not find a valid amplitude for PCHIP demonstration.');
end


%% =========================================================================
% OPTIONAL: PRINT A SHORT Q&A SUMMARY TABLE
% =========================================================================

fprintf('\n============================================================\n');
fprintf('Q&A SUMMARY\n');
fprintf('============================================================\n');

fprintf('\nDuffing model:\n');
fprintf('  Alpha comes from least-squares fit of f_res^2 vs A^2.\n');
fprintf('  R2 quantifies adequacy of the cubic representation.\n');

fprintf('\nBFS:\n');
fprintf('  Uses ONLY first and last valid resonance points.\n');
fprintf('  No Duffing fit enters the BFS calculation.\n');

fprintf('\nPeak identification:\n');
fprintf('  Experimental sweep spacing = %.4f Hz.\n',...
    config.steppedsine.freq.Step);
fprintf('  Resonance is refined using PCHIP interpolation.\n');
fprintf('  PCHIP refines peak location but adds no new experimental information.\n');

fprintf('\nFigures saved to:\n  %s\n\n',qaDir);
%% ==== END QA_BackupPlots.m ====


%% ==== BEGIN PlotBFSvsAlphaComparison.m (integrated) ====
%% PlotBFSvsAlphaComparison.m
% ============================================================
% Side-by-side comparison:
% BFS endpoint metric vs alpha least-squares fit
% for representative temperatures: -30, 20, 70, 120 C
%
% Run AFTER process2.m
% ============================================================

clearvars -except proc_all config

sysLabel   = 'Jointed';
rep_temps  = [-30, 20, 70, 120];

% Recover temperature vector
if isfield(proc_all.(sysLabel).thermo.duffing,'T_vec')
    T_vec = proc_all.(sysLabel).thermo.duffing.T_vec(:);
else
    T_vec = config.steppedsine.temperature.Vector(:);
end

% Output folder
outDir = fullfile('exports','figures','QA');
if ~exist(outDir,'dir')
    mkdir(outDir);
end

% ------------------------------------------------------------
% OPTIONAL: summary figure containing all 4 temperatures
% 4 rows x 2 columns:
% left = BFS, right = alpha
% ------------------------------------------------------------
figSummary = figure('Name','BFS vs Alpha Comparison - Representative Temperatures', ...
    'Color','w', ...
    'Position',[100 50 1200 1400]);

tiledlayout(4,2,'Padding','compact','TileSpacing','compact');

% ------------------------------------------------------------
% Loop over representative temperatures
% ------------------------------------------------------------
for iT = 1:numel(rep_temps)

    T_demo = rep_temps(iT);

    % Find closest available temperature index
    [~, Z] = min(abs(T_vec - T_demo));
    T_actual = T_vec(Z);

    % Backbone data at this temperature
    bb = proc_all.(sysLabel).steppedsine.backbone{Z};

    A_um  = bb.A_pts(:) * 1e6;   % displacement amplitude [um]
    f_res = bb.f_pts(:);         % resonance frequency [Hz]

    valid = ~isnan(A_um) & ~isnan(f_res) & A_um > 0 & f_res > 0;
    A_um  = A_um(valid);
    f_res = f_res(valid);

    if numel(A_um) < 2
        warning('Not enough valid backbone points at %.0f°C. Skipping.', T_actual);
        continue
    end

    % --------------------------------------------------------
    % BFS quantities
    % --------------------------------------------------------
    A_low  = A_um(1);
    A_high = A_um(end);
    f_low  = f_res(1);
    f_high = f_res(end);
    BFS    = (f_high - f_low) / f_low * 100;

    % --------------------------------------------------------
    % Duffing quantities
    % --------------------------------------------------------
    x = A_um.^2;    % A^2 [um^2]
    y = f_res.^2;   % f_res^2 [Hz^2]

    alpha = proc_all.(sysLabel).thermo.duffing.alpha(Z);
    f0_sq = proc_all.(sysLabel).thermo.duffing.f0sq(Z);

    if isfield(proc_all.(sysLabel).thermo.duffing,'R2')
        R2 = proc_all.(sysLabel).thermo.duffing.R2(Z);
    else
        R2 = NaN;
    end

    x_fit = linspace(0, max(x)*1.05, 200);
    y_fit = f0_sq + alpha*x_fit;

    % ========================================================
    % INDIVIDUAL FIGURE FOR THIS TEMPERATURE
    % ========================================================
    fig = figure('Name',sprintf('BFS vs Alpha - %s - %.0fC',sysLabel,T_actual), ...
        'Color','w', ...
        'Position',[100 100 1200 500]);

    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    % ---------------- LEFT: BFS ----------------
    nexttile;
    hold on; box on; grid on;

    plot(A_um, f_res, 'o-', ...
        'LineWidth', 1.6, ...
        'MarkerSize', 7, ...
        'DisplayName', 'Measured resonance points');

    scatter(A_low, f_low, 90, 'filled', ...
        'DisplayName', 'Lowest-amplitude point');

    scatter(A_high, f_high, 90, 'filled', ...
        'DisplayName', 'Highest-amplitude point');

    text(A_low, f_low, sprintf('  f_{low} = %.3f Hz', f_low), ...
        'FontSize', 11, 'VerticalAlignment','bottom');

    text(A_high, f_high, sprintf('  f_{high} = %.3f Hz', f_high), ...
        'FontSize', 11, 'VerticalAlignment','top');

    xlabel('Response displacement amplitude (\mum)', 'FontSize', 12);
    ylabel('Resonance frequency (Hz)', 'FontSize', 12);

    title(sprintf('BFS endpoint metric\nBFS = %.3f%%', BFS), ...
        'FontSize', 13, 'FontWeight','bold');

    legend('Location','best');

    % ---------------- RIGHT: alpha ----------------
    nexttile;
    hold on; box on; grid on;

    plot(x, y, 'o', ...
        'MarkerSize', 7, ...
        'LineWidth', 1.4, ...
        'DisplayName', 'Measured backbone points');

    plot(x_fit, y_fit, '-', ...
        'LineWidth', 2.0, ...
        'DisplayName', 'Least-squares fit');

    scatter(0, f0_sq, 80, 'filled', ...
        'DisplayName', 'Intercept = f_0^2');

    xlabel('A^2 (\mum^2)', 'FontSize', 12);
    ylabel('f_{res}^2 (Hz^2)', 'FontSize', 12);

    title(sprintf('Cubic fit used for \\alpha\n\\alpha = %.3e Hz^2/\\mum^2', alpha), ...
        'FontSize', 13, 'FontWeight','bold');

    x_mid = 0.55 * max(x_fit);
    y_mid = f0_sq + alpha*x_mid;

    text(x_mid, y_mid, '  slope = \alpha', ...
        'FontSize', 11, 'VerticalAlignment','bottom');

    text(0.02*max(x_fit), f0_sq, 'intercept = f_0^2', ...
        'FontSize', 11, 'VerticalAlignment','bottom');

    legend('Location','best');

    sgtitle(sprintf('Why BFS and \\alpha summarise the backbone differently — %s, %.0f°C', ...
        sysLabel, T_actual), ...
        'FontSize', 15, 'FontWeight', 'bold');

    % Save individual figure
    outName = sprintf('QA_BFS_vs_alpha_comparison_%s_T%+03.0fC.png', sysLabel, T_actual);
    exportgraphics(fig, fullfile(outDir, outName), 'Resolution', 300);

    fprintf('\nSaved figure: %s\n', fullfile(outDir, outName));
    fprintf('T = %.0f C | BFS = %.3f %% | alpha = %.3e Hz^2/um^2 | R2 = %.3f\n', ...
        T_actual, BFS, alpha, R2);

    % ========================================================
    % ADD SAME CONTENT TO SUMMARY FIGURE
    % ========================================================
    figure(figSummary);

    % ---- Summary tile: BFS ----
    nexttile;
    hold on; box on; grid on;

    plot(A_um, f_res, 'o-', ...
        'LineWidth', 1.2, ...
        'MarkerSize', 5);

    scatter(A_low, f_low, 60, 'filled');
    scatter(A_high, f_high, 60, 'filled');

    xlabel('A (\mum)', 'FontSize', 10);
    ylabel('f_{res} (Hz)', 'FontSize', 10);

    title(sprintf('%.0f°C: BFS = %.3f%%', T_actual, BFS), ...
        'FontSize', 11, 'FontWeight','bold');

    % ---- Summary tile: alpha ----
    nexttile;
    hold on; box on; grid on;

    plot(x, y, 'o', ...
        'MarkerSize', 5, ...
        'LineWidth', 1.1);

    plot(x_fit, y_fit, '-', ...
        'LineWidth', 1.6);

    scatter(0, f0_sq, 50, 'filled');

    xlabel('A^2 (\mum^2)', 'FontSize', 10);
    ylabel('f_{res}^2 (Hz^2)', 'FontSize', 10);

    title(sprintf('%.0f°C: \\alpha = %.2e, R^2 = %.2f', T_actual, alpha, R2), ...
        'FontSize', 11, 'FontWeight','bold');
end

% ------------------------------------------------------------
% Save summary figure
% ------------------------------------------------------------
figure(figSummary);
sgtitle(sprintf('BFS endpoint metric vs Duffing \\alpha fit — %s', sysLabel), ...
    'FontSize', 16, 'FontWeight', 'bold');

exportgraphics(figSummary, ...
    fullfile(outDir, 'QA_BFS_vs_alpha_summary_4temps.png'), ...
    'Resolution', 300);

fprintf('\nSaved summary figure: %s\n', ...
    fullfile(outDir, 'QA_BFS_vs_alpha_summary_4temps.png'));
%% ==== END PlotBFSvsAlphaComparison.m ====

%% ===== REPRODUCIBLE NUMERICAL VALIDATION EXPORTS =====
% Do not modify any computed fitting values: export from proc_all only.
validationDir = fullfile('exports', 'validation');
if ~exist(validationDir, 'dir'), mkdir(validationDir); end
validationTimestamp = datestr(now, 'yyyy-mm-ddTHH-MM-SS');
validationSystems = {'Monolithic','Jointed'};
validationTemperature_C = (-30:10:120).';
coeffSystem = strings(0,1); coeffT = zeros(0,1);
coeffAlpha = zeros(0,1); coeffF0sq = zeros(0,1);
coeffR2 = zeros(0,1); coeffDegree = zeros(0,1);
polySystem = strings(0,1); polyQuantity = strings(0,1);
polyDegree = zeros(0,1); polyR2 = zeros(0,1); polyAdjustedR2 = zeros(0,1);
backboneSystem = strings(0,1); backboneT = zeros(0,1);
backboneAmplitude_m = zeros(0,1); backboneFrequency_Hz = zeros(0,1);
backboneIndex = zeros(0,1);
for validationSystemIndex = 1:numel(validationSystems)
    validationSystem = validationSystems{validationSystemIndex};
    validationDuffing = proc_all.(validationSystem).thermo.duffing;
    nTemperatures = numel(validationDuffing.T_vec);
    if nTemperatures ~= 16
        error('ResearchFinal:TemperatureCount', ...
            '%s: expected 16 fitted temperatures; found %d.',validationSystem,nTemperatures);
    end
    thisT = validationDuffing.T_vec(:);
    thisAlpha = validationDuffing.alpha(:);
    thisF0sq = validationDuffing.f0sq(:);
    thisR2 = validationDuffing.R2(:);
    if numel(thisAlpha) ~= nTemperatures || numel(thisF0sq) ~= nTemperatures || numel(thisR2) ~= nTemperatures
        error('ResearchFinal:CoefficientShape','Unexpected Duffing array sizes for %s.',validationSystem);
    end
    coeffSystem = [coeffSystem; repmat(string(validationSystem), nTemperatures, 1)]; %#ok<AGROW>
    coeffT = [coeffT; thisT]; %#ok<AGROW>
    coeffAlpha = [coeffAlpha; thisAlpha]; %#ok<AGROW>
    coeffF0sq = [coeffF0sq; thisF0sq]; %#ok<AGROW>
    coeffR2 = [coeffR2; thisR2]; %#ok<AGROW>
    coeffDegree = [coeffDegree; repmat(validationDuffing.best_deg,nTemperatures,1)]; %#ok<AGROW>
    polyFields = {'alpha_poly','f0sq_poly'};
    polyLabels = {'Duffing_alpha','f0_squared'};
    for quantityIndex = 1:2
        candidates = validationDuffing.(polyFields{quantityIndex});
        for degreeIndex = 1:numel(candidates)
            polySystem(end+1,1) = string(validationSystem); %#ok<SAGROW>
            polyQuantity(end+1,1) = string(polyLabels{quantityIndex}); %#ok<SAGROW>
            polyDegree(end+1,1) = degreeIndex; %#ok<SAGROW>
            polyR2(end+1,1) = candidates(degreeIndex).R2; %#ok<SAGROW>
            polyAdjustedR2(end+1,1) = candidates(degreeIndex).adjR2; %#ok<SAGROW>
        end
    end
    validationBackbones = proc_all.(validationSystem).steppedsine.backbone;
    for temperatureIndex = 1:numel(validationBackbones)
        b = validationBackbones{temperatureIndex};
        frequencyValues = b.f_pts(:);
        amplitudeValues = b.A_pts(:);
        if numel(frequencyValues) ~= numel(amplitudeValues)
            error('ResearchFinal:BackboneShape', ...
                'Mismatched backbone arrays for %s, temperature index %d.',validationSystem,temperatureIndex);
        end
        nPoints = numel(frequencyValues);
        backboneSystem = [backboneSystem; repmat(string(validationSystem),nPoints,1)]; %#ok<AGROW>
        backboneT = [backboneT; repmat(thisT(temperatureIndex),nPoints,1)]; %#ok<AGROW>
        backboneAmplitude_m = [backboneAmplitude_m; amplitudeValues]; %#ok<AGROW>
        backboneFrequency_Hz = [backboneFrequency_Hz; frequencyValues]; %#ok<AGROW>
        backboneIndex = [backboneIndex; (1:nPoints).']; %#ok<AGROW>
    end
end
coefficientTable = table(coeffSystem,coeffT,coeffAlpha,coeffF0sq,coeffR2,coeffDegree, ...
    'VariableNames',{'System','Temperature_C','Duffing_alpha_Hz2_per_um2', ...
    'f0_squared_Hz2','DuffingFit_R2','SelectedPolynomialDegree'});
polynomialTable = table(polySystem,polyQuantity,polyDegree,polyR2,polyAdjustedR2, ...
    'VariableNames',{'System','Quantity','Degree','R2','AdjustedR2'});
backboneTable = table(backboneSystem,backboneT,backboneIndex,backboneAmplitude_m,backboneFrequency_Hz, ...
    'VariableNames',{'System','Temperature_C','AmplitudeIndex','DisplacementAmplitude_m','ResonanceFrequency_Hz'});
writetable(coefficientTable,fullfile(validationDir,'duffing_temperature_coefficients.csv'));
writetable(polynomialTable,fullfile(validationDir,'temperature_polynomial_candidates.csv'));
writetable(backboneTable,fullfile(validationDir,'resonance_backbone_points.csv'));
%% ===== OPTIONAL HIGHER-ORDER BACKBONE MODEL COMPARISON =====
% This is a secondary, per-temperature comparison. The existing first-order
% Duffing coefficients, thermal polynomials, and baseline plots are unchanged.
% Extended model: f_res^2 = f0sq + alpha*A_um^2 + beta4*A_um^4.
% Model selection uses out-of-sample LOOCV frequency RMSE, with an AICc check.
comparisonOutputDir = fullfile(validationDir, 'higher_order_backbone');
if ~isfolder(comparisonOutputDir), mkdir(comparisonOutputDir); end
[extendedModelTable, extendedPredictionTable] = ...
    compareHigherOrderBackbones(backboneTable,comparisonOutputDir);
writetable(extendedModelTable, fullfile(comparisonOutputDir,'backbone_model_comparison.csv'));
writetable(extendedPredictionTable, fullfile(comparisonOutputDir,'backbone_model_predictions.csv'));

%% ===== EXTENDED MODEL IDENTIFIABILITY AND TURNING-POINT VALIDATION =====
% Diagnostic only: does NOT change the original Duffing parameters or
% the exploratory model-preference decisions above.
[uncertaintyTable, influenceTable] = validateExtendedBackboneModels(backboneTable, extendedModelTable);
writetable(uncertaintyTable,fullfile(comparisonOutputDir,'extended_parameter_validation.csv'));
writetable(influenceTable,fullfile(comparisonOutputDir,'extended_leave_one_out_influence.csv'));
plotColdJointedValidation(backboneTable,extendedModelTable,uncertaintyTable,comparisonOutputDir);

%% ===== PEAK-RESOLUTION SENSITIVITY OF EXTENDED MODEL =====
% Diagnostic only. Re-evaluates raw-grid / PCHIP / quadratic peak estimators
% while holding displacement amplitudes fixed. This does not change the
% original Duffing coefficients or model screening.
peakSensitivityDir = fullfile(comparisonOutputDir,'peak_resolution');
if ~isfolder(peakSensitivityDir), mkdir(peakSensitivityDir); end
[peakSensitivityPoints, peakSensitivityComparison] = ...
    validatePeakResolutionHigherOrder(proc_all,peakSensitivityDir);

fprintf('Extended-model identifiability outputs written to %s\n',comparisonOutputDir);

fprintf('Higher-order backbone comparison written to %s\\n',comparisonOutputDir);

% MATLAB numerical snapshot only, avoiding all transient large work arrays.
save(fullfile(validationDir,'research_results_compact.mat'), 'proc_all','config','-v7.3');
manifestFile = fopen(fullfile(validationDir,'RUN_MANIFEST.txt'),'w');
if manifestFile < 0, error('ResearchFinal:CannotWriteManifest','Could not open validation manifest.'); end
fprintf(manifestFile,'Run timestamp: %s\n',validationTimestamp);
fprintf(manifestFile,'Pipeline: SUTMNT_experimental_pipeline_release.m\n');
fprintf(manifestFile,'MATLAB: %s\n',version);
fprintf(manifestFile,'Systems: Monolithic, Jointed\n');
fprintf(manifestFile,'Duffing observations: %d\n',height(coefficientTable));
fprintf(manifestFile,'Backbone points: %d\n',height(backboneTable));
fprintf(manifestFile,'Selected polynomial degree is based on alpha adjusted R2 and applied to f0 squared.\n');
fprintf(manifestFile,'Important: stiffness-related and energy-proxy quantities are diagnostic; physical energy dissipation is not calibrated.\n');
fprintf(manifestFile,'Peak estimation: PCHIP baseline; quadratic alternative was a sensitivity test only.\n');
fprintf(manifestFile,'Interpret Jointed -30/-20/-10 C cubic coefficients cautiously because those fits are weakly identified.\n');
fprintf(manifestFile,'Higher-order backbone comparison is exploratory; primary Duffing coefficients unchanged.\n');
fclose(manifestFile);
fprintf('Validation exports written to %s\n',validationDir);

%% ===== JOURNAL FIGURES AND TABLES =====
% The journal output generator reads ONLY the validated exports written above.
% It does not modify the scientific processing, fitted parameters, or exports.
journalGenerator = fullfile(pwd,'generateJournalOutputs.m'); % final publication generator
if ~isfile(journalGenerator)
    error('ResearchFinal:MissingJournalGenerator', ...
        ['generateJournalOutputs.m was not found in the project folder. ' ...
         'Keep it beside SUTMNT_experimental_pipeline.m.']);
end

fprintf('\nGenerating journal figures and tables from validated exports...\n');
run(journalGenerator);

% Verify every numerically generated main-paper figure and table exists.
% Figure 01 (specimen/test-setup artwork) is intentionally external/manual.
journalDir = fullfile('exports','journal');
journalExpected = { ...
    'Figure_02_impact_FRF.pdf', ...
    'Figure_03_backbone_evolution.pdf', ...
    'Figure_04_20C_comparison.pdf', ...
    'Figure_05_jointed_cold_backbones.pdf', ...
    'Figure_06_BFS_and_alpha.pdf', ...
    'Figure_07_standard_vs_extended.pdf', ...
    'Figure_08_model_validation.pdf', ...
    'Figure_09_cold_jointed_turning_point_context.pdf', ...
    'Table_01_baseline.csv', ...
    'Table_02_temperature_results.csv', ...
    'Table_03_extended_validation.csv', ...
    'manuscript_numbers.txt'};

missingJournalOutputs = strings(0,1);
for iJournal = 1:numel(journalExpected)
    thisJournalFile = fullfile(journalDir,journalExpected{iJournal});
    if ~isfile(thisJournalFile)
        missingJournalOutputs(end+1,1) = string(thisJournalFile); %#ok<AGROW>
    end
end

if ~isempty(missingJournalOutputs)
    error('ResearchFinal:MissingJournalOutputs', ...
        'Journal output generation completed but %d expected files are missing:\n%s', ...
        numel(missingJournalOutputs), strjoin(missingJournalOutputs,newline));
end

fprintf('Journal output verification passed: all generated paper figures/tables are present.\n');
fprintf('Note: Figure_01_test_setup.pdf is manual artwork and is not generated by the numerical pipeline.\n');

diary off;

%% ===== LOCAL HELPER FUNCTIONS =====

% ---- Higher-order backbone analysis: independent of primary Duffing fits ----
function [detail,summary] = validatePeakResolutionHigherOrder(proc_all,outputDir)
%VALIDATEPEAKRESOLUTIONHIGHERORDER Compare raw-grid, PCHIP and quadratic peak localisation.
%   This is a numerical sensitivity check only; it is not repeatability testing.
% Compare peak estimators for the extended backbone. Keep original A fixed to
% isolate peak-frequency extraction. No changes to primary parameters.
systems = {'Monolithic','Jointed'};
detail = table(); summary = table();
methods = {'RawGrid','PCHIP','Quadratic'};
for si = 1:numel(systems)
    sys = systems{si};
    saved = proc_all.(sys).steppedsine;
    temps = proc_all.(sys).thermo.duffing.T_vec(:);
    for ti = 1:numel(temps)
        T = temps(ti);
        bb = saved.backbone{ti};
        amp = saved.amp{ti};
        n = numel(bb.A_pts);
        A = bb.A_pts(:) * 1e6; % um
        F = NaN(n,3);
        for ai = 1:min(n,numel(amp))
            if ~isfield(amp(ai),'Fit') || ~isfield(amp(ai).Fit,'fRaw') || ...
                    ~isfield(amp(ai).Fit,'yRaw'), continue; end
            fx = amp(ai).Fit.fRaw(:);
            yy = amp(ai).Fit.yRaw(:);
            ok = isfinite(fx)&isfinite(yy);
            fx=fx(ok); yy=yy(ok);
            if numel(fx)<3, continue; end
            [fx,ix] = sort(fx); yy=yy(ix);
            [fx,ix] = unique(fx,'stable'); yy=yy(ix);
            if numel(fx)<3,continue;end
            [~,pk] = max(yy);
            F(ai,1)=fx(pk);
            % Use the actual pipeline result for PCHIP: avoids discrepancies
            % caused by reimplementing its dense interpolation grid.
            F(ai,2)=bb.f_pts(ai);
            if pk>1 && pk<numel(fx)
                xx=fx(pk-1:pk+1); zz=yy(pk-1:pk+1);
                pp=polyfit(xx-xx(2),zz,2);
                if pp(1)<0
                    shift=-pp(2)/(2*pp(1));
                    if isfinite(shift)&& abs(shift)<=max(diff(xx))
                        F(ai,3)=xx(2)+shift;
                    end
                end
            end
            if ~isfinite(F(ai,3)),F(ai,3)=F(ai,1);end
        end
        for mi=1:3
            [fit,ok] = fitPeakMethodHigherOrder(A,F(:,mi));
            if ~ok,continue;end
            row = table(string(sys),T,string(methods{mi}),fit.n, ...
                fit.rmseStd,fit.rmseExt,fit.cvStd,fit.cvExt,fit.aiccStd,fit.aiccExt, ...
                fit.turnMM,fit.turnInside,fit.alpha,fit.beta4, ...
                'VariableNames',{'System','Temperature_C','Method','N', ...
                'StandardRMSE_Hz','ExtendedRMSE_Hz','StandardLOOCV_Hz', ...
                'ExtendedLOOCV_Hz','StandardAICc','ExtendedAICc', ...
                'TurningPoint_mm','TurningPointInRange','ExtendedAlpha_Hz2_um2', ...
                'ExtendedBeta4_Hz2_um4'});
            summary=[summary;row]; %#ok<AGROW>
            good = isfinite(A)&isfinite(F(:,mi))&A>0&F(:,mi)>0;
            aidx=find(good);
            rowDetail=table(repmat(string(sys),numel(aidx),1),repmat(T,numel(aidx),1), ...
                repmat(string(methods{mi}),numel(aidx),1),aidx, ...
                A(good)/1000,F(good,mi), ...
                'VariableNames',{'System','Temperature_C','Method', ...
                'AmplitudeIndex','Amplitude_mm','PeakFrequency_Hz'});
            detail=[detail;rowDetail]; %#ok<AGROW>
        end
    end
end
writetable(detail,fullfile(outputDir,'peak_estimator_backbone_points.csv'));
writetable(summary,fullfile(outputDir,'peak_estimator_extended_comparison.csv'));
% Dedicated turning-point overview. Plot only model-derived values.
fig=figure('Visible','off','Color','w','Position',[100 100 1100 350]);
tiledlayout(fig,1,3,'TileSpacing','compact');
for i=1:3
    T=[-30,-20,-10];
    ax=nexttile; hold(ax,'on');
    rr=summary.System=="Jointed" & summary.Temperature_C==T(i);
    ss=summary(rr,:);
    bar(ax,categorical(ss.Method,methods),ss.TurningPoint_mm,'FaceColor',[.35 .45 .65]);
    yline(ax,0,'k-');
    ylabel(ax,'Turning-point amplitude (mm)');
    title(ax,sprintf('Jointed, %d °C',T(i)));
    grid(ax,'on');
end
try
    exportgraphics(fig,fullfile(outputDir,'peak_method_turning_points.png'),'Resolution',300);
    exportgraphics(fig,fullfile(outputDir,'peak_method_turning_points.pdf'),'ContentType','vector');
catch ME
    warning('PeakResolution:Export','Plot export failed: %s',ME.message);
end
close(fig);
fprintf('Peak-resolution sensitivity: %d fits exported to %s\n',height(summary),outputDir);
end

function [fit,ok] = fitPeakMethodHigherOrder(A,f)
%FITPEAKMETHODHIGHERORDER Fit the empirical A^2 + A^4 backbone model for one peak method.
fit=struct('n',NaN,'rmseStd',NaN,'rmseExt',NaN, ...
    'cvStd',NaN,'cvExt',NaN,'aiccStd',NaN,'aiccExt',NaN, ...
    'turnMM',NaN,'turnInside',false,'alpha',NaN,'beta4',NaN);
good=isfinite(A)&isfinite(f)&A>0&f>0;
A=A(good);f=f(good); n=numel(f);ok=false;
if n<6||numel(unique(A))<4,return;end
scale=max(A);z=(A/scale).^2;
X2=[ones(n,1),z]; X4=[X2,z.^2];
if rank(X4)<3,return;end
y=f.^2;b2=X2\y;b4=X4\y;
p2=sqrt(max(X2*b2,0));p4=sqrt(max(X4*b4,0));
fit.n=n;fit.rmseStd=sqrt(mean((f-p2).^2));
fit.rmseExt=sqrt(mean((f-p4).^2));
% Explicit leave-one-out refits, evaluated in the frequency domain.
e2=NaN(n,1);e4=NaN(n,1);
for k=1:n
    tr=true(n,1);tr(k)=false;
    if rank(X4(tr,:))<3,continue;end
    b2k=X2(tr,:)\y(tr); b4k=X4(tr,:)\y(tr);
    e2(k)=f(k)-sqrt(max(X2(k,:)*b2k,0));
    e4(k)=f(k)-sqrt(max(X4(k,:)*b4k,0));
end
fit.cvStd=sqrt(mean(e2.^2,'omitnan'));
fit.cvExt=sqrt(mean(e4.^2,'omitnan'));
% Comparable AICc computed on squared-frequency residuals.
sse2=max(sum((y-X2*b2).^2),realmin);
sse4=max(sum((y-X4*b4).^2),realmin);
fit.aiccStd=n*log(sse2/n)+2*2+2*2*3/(n-2-1);
fit.aiccExt=n*log(sse4/n)+2*3+2*3*4/(n-3-1);
fit.alpha=b4(2)/scale^2; fit.beta4=b4(3)/scale^4;
if b4(3)>0 && b4(2)<0
    az=sqrt(-b4(2)/(2*b4(3)));
    fit.turnMM=(scale*az)/1000;
    fit.turnInside=fit.turnMM>=min(A)/1000&&fit.turnMM<=max(A)/1000;
end
ok=true;
end

function [summary, predictions] = compareHigherOrderBackbones(points,outputDir)
%COMPAREHIGHERORDERBACKBONES Compare standard A^2 and extended A^2+A^4 models by temperature.
%   Selection uses LOOCV frequency RMSE together with AICc screening.
% Compare standard vs extended Duffing using *identical* valid measurements.
% Scaled coordinate z=(A_um / max(A_um))^2 guards A^4 conditioning.
% Coefficients are returned in Hz^2/um^2 and Hz^2/um^4.
systemNames = ["Monolithic","Jointed"];
tempCandidates = unique(points.Temperature_C(:));
summary = table();
predictions = table();
for si = 1:numel(systemNames)
    systemName = systemNames(si);
    for ti = 1:numel(tempCandidates)
        T = tempCandidates(ti);
        rows = string(points.System)==systemName & points.Temperature_C==T;
        A = points.DisplacementAmplitude_m(rows) * 1e6; % micrometres
        f = points.ResonanceFrequency_Hz(rows);
        pointIndex = points.AmplitudeIndex(rows);
        valid = isfinite(A) & isfinite(f) & A>0 & f>0;
        A = A(valid); f = f(valid); pointIndex = pointIndex(valid);
        n = numel(f);
        if n < 6 || numel(unique(A)) < 3
            warning('HigherOrder:InsufficientPoints', ...
                '%s %g C: fewer than six usable observations.',systemName,T);
            continue;
        end
        scaleUm = max(A);
        z = (A/scaleUm).^2;
        y = f.^2;
        X2 = [ones(n,1), z];
        X4 = [X2, z.^2];
        if rank(X4)<3
            warning('HigherOrder:RankDeficient','%s %g C: A^4 model not identifiable.',systemName,T);
            continue;
        end
        b2 = X2\y; b4 = X4\y;
        pred2 = sqrt(max(X2*b2,0));
        pred4 = sqrt(max(X4*b4,0));
        rmse2 = sqrt(mean((f-pred2).^2));
        rmse4 = sqrt(mean((f-pred4).^2));
        r2_2 = backboneR2(y,X2*b2);
        r2_4 = backboneR2(y,X4*b4);
        [cv2, cv4] = backboneLOOCV(z, f);
        aicc2 = backboneAICc(y,X2*b2,2);
        aicc4 = backboneAICc(y,X4*b4,3);
        % AICc support AND at least 5% improvement in leave-one-out RMSE.
        % This is a predeclared exploratory screening rule, not a significance test.
        supportsExtended = isfinite(cv2) && isfinite(cv4) && cv4 <= 0.95*cv2 && ...
            aicc4 <= aicc2-2;
        if supportsExtended
            selected = "Extended";
        else
            selected = "Standard";
        end
        relImprovementPercent = 100*(cv2-cv4)/max(cv2,eps);
        summaryRow = table(systemName,T,n,b2(1),b2(2)/scaleUm^2, ...
            b4(1),b4(2)/scaleUm^2,b4(3)/scaleUm^4, ...
            r2_2,r2_4,rmse2,rmse4,cv2,cv4,aicc2,aicc4, ...
            relImprovementPercent,selected, ...
            'VariableNames',{'System','Temperature_C','N','Standard_f0sq_Hz2', ...
            'Standard_alpha_Hz2_per_um2','Extended_f0sq_Hz2', ...
            'Extended_alpha_Hz2_per_um2','Extended_beta4_Hz2_per_um4', ...
            'Standard_R2_fSquared','Extended_R2_fSquared','Standard_TrainRMSE_Hz', ...
            'Extended_TrainRMSE_Hz','Standard_LOOCV_RMSE_Hz', ...
            'Extended_LOOCV_RMSE_Hz','Standard_AICc','Extended_AICc', ...
            'LOOCV_Improvement_pct','ScreenedPreference'});
        summary = [summary;summaryRow]; %#ok<AGROW>
        predictionRows = table(repmat(systemName,n,1),repmat(T,n,1),pointIndex,A/1000, ...
            f,pred2,pred4,'VariableNames',{'System','Temperature_C', ...
            'AmplitudeIndex','Amplitude_mm','Observed_Hz','StandardPredicted_Hz', ...
            'ExtendedPredicted_Hz'});
        predictions = [predictions;predictionRows]; %#ok<AGROW>
    end
end
if isempty(summary)
    warning('HigherOrder:NoModels','No higher-order backbone models could be fitted.');
    return;
end
% Diagnostic plot: cold jointed specimen plus 20 C. No fabricated data.
fig = figure('Visible','off','Color','w','Units','inches','Position',[1 1 10 4]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
for ix=1:2
    T = [-30,20]; T=T(ix);
    ax=nexttile(tl);hold(ax,'on');
    row=predictions(string(predictions.System)=="Jointed" & predictions.Temperature_C==T,:);
    if isempty(row), title(ax,sprintf('Jointed, %d C: no fit',T));continue;end
    [x,idx]=sort(row.Amplitude_mm);
    y=row.Observed_Hz(idx);
    plot(ax,x,y,'ko','MarkerSize',5,'DisplayName','Measured');
    plot(ax,x,row.StandardPredicted_Hz(idx),'b-','LineWidth',1.4,'DisplayName','Standard Duffing');
    plot(ax,x,row.ExtendedPredicted_Hz(idx),'r--','LineWidth',1.4,'DisplayName','Extended A^4');
    grid(ax,'on');xlabel(ax,'Displacement amplitude (mm)');
    ylabel(ax,'Resonance frequency (Hz)');title(ax,sprintf('Jointed at %+d C',T));
    legend(ax,'Location','best','Box','off');
end
try
    exportgraphics(fig,fullfile(outputDir,'jointed_backbone_model_comparison.png'),'Resolution',300);
    exportgraphics(fig,fullfile(outputDir,'jointed_backbone_model_comparison.pdf'),'ContentType','vector');
catch ME
    warning('HigherOrder:PlotExport','Higher-order comparison plot could not be exported: %s',ME.message);
end
close(fig);
fig=figure('Visible','off','Color','w','Units','inches','Position',[1 1 9 4.5]);
t= tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
for si=1:2
    ax=nexttile(t);hold(ax,'on');
    r=summary(string(summary.System)==systemNames(si),:);
    plot(ax,r.Temperature_C,r.Standard_LOOCV_RMSE_Hz,'o-','DisplayName','Standard');
    plot(ax,r.Temperature_C,r.Extended_LOOCV_RMSE_Hz,'s--','DisplayName','Extended A^4');
    xlabel(ax,'Temperature (C)');ylabel(ax,'LOOCV RMSE (Hz)');
    title(ax,char(systemNames(si)));grid(ax,'on');legend(ax,'Location','best','Box','off');
end
try
    exportgraphics(fig,fullfile(outputDir,'backbone_LOOCV_by_temperature.png'),'Resolution',300);
    exportgraphics(fig,fullfile(outputDir,'backbone_LOOCV_by_temperature.pdf'),'ContentType','vector');
catch ME
    warning('HigherOrder:PlotExport','LOOCV plot could not be exported: %s',ME.message);
end
close(fig);
end

function [summary, influence] = validateExtendedBackboneModels(points,modelTable)
%VALIDATEEXTENDEDBACKBONEMODELS Quantify conditional fit sensitivity of the extended model.
%   Includes residual bootstrap, parameter correlation and leave-one-out influence diagnostics.
% Confidence intervals describe conditional least-squares fitting uncertainty.
% They do NOT include calibrated sensor, repeatability, or peak-estimator error.
% Resampling is residual bootstrap at fixed measured amplitudes, B=500.
rng(17031,'twister');
B=500;
summary=table(); influence=table();
for ri=1:height(modelTable)
    row=modelTable(ri,:);
    sys=string(row.System); T=row.Temperature_C;
    mask=string(points.System)==sys & points.Temperature_C==T;
    A=points.DisplacementAmplitude_m(mask)*1e6;
    f=points.ResonanceFrequency_Hz(mask);
    valid=isfinite(A)&isfinite(f)&A>0&f>0;
    A=A(valid);f=f(valid);
    n=numel(f);
    if n<8 || numel(unique(A))<4,continue;end
    scale=max(A);z=(A./scale).^2;
    X=[ones(n,1),z,z.^2]; y=f.^2;
    if rank(X)<3,continue;end
    b=X\y; yh=X*b; residual=y-yh;
    sigma2=sum(residual.^2)/(n-3);
    covariance=sigma2*((X'*X)\eye(3));
    covariance=(covariance+covariance')/2;
    se=sqrt(max(diag(covariance),0));
    corrAlphaBeta=covariance(2,3)/max(se(2)*se(3),eps);
    alpha=b(2)/scale^2; beta4=b(3)/scale^4;
    alphaSE=se(2)/scale^2;betaSE=se(3)/scale^4;
    turnUm=NaN; turnWithin=false; turnType="None";
    if beta4~=0 && -alpha/(2*beta4)>0
        turnUm=sqrt(-alpha/(2*beta4));
        turnWithin=(turnUm>=min(A) && turnUm<=max(A));
        if beta4>0, turnType="Minimum"; else, turnType="Maximum";end
    end
    bootAlpha=NaN(B,1);bootBeta=NaN(B,1);bootTurn=NaN(B,1);
    bootInside=false(B,1);
    for k=1:B
        yb=yh+residual(randi(n,n,1));
        bb=X\yb;
        bootAlpha(k)=bb(2)/scale^2;
        bootBeta(k)=bb(3)/scale^4;
        if bootBeta(k)~=0 && -bootAlpha(k)/(2*bootBeta(k))>0
            tb=sqrt(-bootAlpha(k)/(2*bootBeta(k)));
            if isfinite(tb)
                bootTurn(k)=tb;
                bootInside(k)=tb>=min(A) && tb<=max(A);
            end
        end
    end
    ciAlpha=prctile(bootAlpha,[2.5,97.5]);
    ciBeta=prctile(bootBeta,[2.5,97.5]);
    validTurn=bootTurn(isfinite(bootTurn));
    ciTurn=[NaN,NaN];
    if numel(validTurn)>=ceil(.5*B)
        ciTurn=prctile(validTurn,[2.5,97.5]);
    end
    fracInside=mean(bootInside);
    % Deterministic, explicitly conservative descriptive flag, not a p-value.
    turnStable=turnWithin && fracInside>=0.8 && ...
        all(isfinite(ciTurn)) && ciTurn(1)>=min(A) && ciTurn(2)<=max(A);
    LOOalpha=NaN(n,1);LOObeta=NaN(n,1);LOOturn=NaN(n,1);
    for k=1:n
        keep=true(n,1);keep(k)=false;
        if rank(X(keep,:))<3,continue;end
        bk=X(keep,:)\y(keep);
        LOOalpha(k)=bk(2)/scale^2;
        LOObeta(k)=bk(3)/scale^4;
        if LOObeta(k)~=0 && -LOOalpha(k)/(2*LOObeta(k))>0
            LOOturn(k)=sqrt(-LOOalpha(k)/(2*LOObeta(k)));
        end
    end
    LOOalphaRange=max(LOOalpha,[],'omitnan')-min(LOOalpha,[],'omitnan');
    LOObetaRange=max(LOObeta,[],'omitnan')-min(LOObeta,[],'omitnan');
    sumrow=table(sys,T,n,min(A)/1000,max(A)/1000,alpha,beta4,alphaSE,betaSE, ...
        ciAlpha(1),ciAlpha(2),ciBeta(1),ciBeta(2),corrAlphaBeta,cond(X), ...
        turnUm/1000,turnWithin,turnType,ciTurn(1)/1000,ciTurn(2)/1000, ...
        fracInside,turnStable,LOOalphaRange,LOObetaRange, ...
        string(row.ScreenedPreference), ...
        'VariableNames',{'System','Temperature_C','N','AmplitudeMin_mm','AmplitudeMax_mm', ...
        'Alpha_Hz2_per_um2','Beta4_Hz2_per_um4','Alpha_SE','Beta4_SE', ...
        'Alpha_BootCI_Lo','Alpha_BootCI_Hi','Beta4_BootCI_Lo','Beta4_BootCI_Hi', ...
        'Alpha_Beta_Correlation','ScaledDesignCondition','TurningPoint_mm', ...
        'TurningPointInRange','TurningType','Turning_BootCI_Lo_mm', ...
        'Turning_BootCI_Hi_mm','BootstrapFractionTurningInRange', ...
        'TurningPointStableScreen','LOO_Alpha_Range','LOO_Beta4_Range','PreviousScreenedPreference'});
    summary=[summary;sumrow]; %#ok<AGROW>
    for k=1:n
        rowk=table(sys,T,k,A(k)/1000,LOOalpha(k),LOObeta(k),LOOturn(k)/1000, ...
            'VariableNames',{'System','Temperature_C','LeftOutIndex','LeftOutAmplitude_mm', ...
            'LOO_Alpha_Hz2_per_um2','LOO_Beta4_Hz2_per_um4','LOO_TurningPoint_mm'});
        influence=[influence;rowk]; %#ok<AGROW>
    end
end
end

function plotColdJointedValidation(points,models,validation,out)
%PLOTCOLDJOINTEDVALIDATION Plot cold jointed fits and model-derived stationary points.
temps=[-30 -20 -10];
fig=figure('Visible','off','Color','w','Units','inches','Position',[1 1 12 4], ...
    'PaperPositionMode','auto');
tl=tiledlayout(fig,1,3,'Padding','compact','TileSpacing','compact');
for ix=1:3
    T=temps(ix);ax=nexttile(tl);hold(ax,'on');
    mask=string(points.System)=="Jointed" & points.Temperature_C==T;
    A=points.DisplacementAmplitude_m(mask)*1e6;
    f=points.ResonanceFrequency_Hz(mask);
    maskM=string(models.System)=="Jointed" & models.Temperature_C==T;
    maskV=string(validation.System)=="Jointed" & validation.Temperature_C==T;
    if ~any(mask) || ~any(maskM),title(ax,sprintf('%d °C: missing',T));continue;end
    [A,order]=sort(A);f=f(order);
    m=models(find(maskM,1),:);
    u=linspace(min(A),max(A),200);
    f2=sqrt(max(m.Standard_f0sq_Hz2+m.Standard_alpha_Hz2_per_um2*u.^2,0));
    f4=sqrt(max(m.Extended_f0sq_Hz2+m.Extended_alpha_Hz2_per_um2*u.^2 + ...
        m.Extended_beta4_Hz2_per_um4*u.^4,0));
    plot(ax,A/1000,f,'ko','MarkerFaceColor','w','MarkerSize',5,'DisplayName','Measured');
    plot(ax,u/1000,f2,'b-','LineWidth',1.4,'DisplayName','Standard');
    plot(ax,u/1000,f4,'r--','LineWidth',1.4,'DisplayName','Extended');
    if any(maskV)
        v=validation(find(maskV,1),:);
        if v.TurningPointInRange && isfinite(v.TurningPoint_mm)
            xline(ax,v.TurningPoint_mm,':','Color',[.5 .5 .5], ...
                'HandleVisibility','off');
        end
    end
    xlabel(ax,'Amplitude (mm)');ylabel(ax,'Frequency (Hz)');
    title(ax,sprintf('Jointed, %d °C',T));grid(ax,'on');
    if ix==1,legend(ax,'Location','best','Box','off');end
end
try
    exportgraphics(fig,fullfile(out,'jointed_low_temperature_validation.png'),'Resolution',300);
    try
        exportgraphics(fig,fullfile(out,'jointed_low_temperature_validation.pdf'),'ContentType','vector');
    catch ME
        warning('ExtendedValidation:PDFExport','PDF failed: %s',ME.message);
    end
catch ME
    warning('ExtendedValidation:FigureExport','Figure export failed: %s',ME.message);
end
close(fig);
end

function r = backboneR2(y,yhat)
%BACKBONER2 Coefficient of determination for backbone predictions.
sst=sum((y-mean(y)).^2);
if sst<=eps, r=NaN;else,r=1-sum((y-yhat).^2)/sst;end
end

function score = backboneAICc(y,yhat,k)
%BACKBONEAICC Small-sample corrected Akaike information criterion.
n=numel(y);
sse=sum((y-yhat).^2);
score = n*log(max(sse/n,realmin)) + 2*k + 2*k*(k+1)/(n-k-1);
end

function [rmse2,rmse4] = backboneLOOCV(z,f)
%BACKBONELOOCV Leave-one-out frequency RMSE for standard and extended backbone models.
n=numel(f);
errors2=NaN(n,1);errors4=NaN(n,1);
for i=1:n
    keep=true(n,1);keep(i)=false;
    X2=[ones(n-1,1),z(keep)];
    X4=[X2,z(keep).^2];
    if rank(X4)<3 || rank(X2)<2,continue;end
    y=f(keep).^2;
    b2=X2\y;b4=X4\y;
    p2=sqrt(max([1,z(i)]*b2,0));
    p4=sqrt(max([1,z(i),z(i).^2]*b4,0));
    errors2(i)=f(i)-p2;errors4(i)=f(i)-p4;
end
rmse2=sqrt(mean(errors2.^2,'omitnan'));
rmse4=sqrt(mean(errors4.^2,'omitnan'));
end

% ---- cycleDeltaFCallback.m ----
function cycleDeltaFCallback(fig,evt)
%CYCLEDELTAFCALLBACK Keyboard callback for interactive diagnostic browsing.

idx   = getappdata(fig,'activeIdx');
temps = getappdata(fig,'temps');
nT    = numel(temps);

switch evt.Key
    case 'rightarrow'
        idx = min(idx+1,nT);
    case 'leftarrow'
        idx = max(idx-1,1);
    otherwise
        return
end

setappdata(fig,'activeIdx',idx);
drawDeltaF(fig);
end

% ---- drawDeltaF.m ----
function drawDeltaF(fig)
%DRAWDELTAF Redraw the active temperature in the Delta-f diagnostic figure.

ax     = getappdata(fig,'ax');
data   = getappdata(fig,'DeltaF');
temps  = getappdata(fig,'temps');
idx    = getappdata(fig,'activeIdx');

cla(ax); hold(ax,'on');

% ---- All temperatures (faint) ----
for i = 1:numel(data)
    plot(ax, data{i}.a, data{i}.J, 'r:', 'LineWidth', 0.8);
    plot(ax, data{i}.a, data{i}.M, 'b:', 'LineWidth', 0.8);
end

% ---- Active temperature (bold) ----
plot(ax, data{idx}.a, data{idx}.J, 'r-', 'LineWidth', 2.5);
plot(ax, data{idx}.a, data{idx}.M, 'b-', 'LineWidth', 2.5);

title(ax, sprintf('Temperature = %d ^\\circC   (← / →)', temps(idx)));
drawnow;
end

% ---- drawFRF.m ----
function drawFRF(figFRF)
%DRAWFRF Redraw stepped-sine FRF diagnostics for the active temperature.

    ax     = getappdata(figFRF,'ax');
    proc   = getappdata(figFRF,'proc');
    temps  = getappdata(figFRF,'temps');
    Na     = getappdata(figFRF,'Na');

    figMain = getappdata(figFRF,'figMain');
    activeZ = getappdata(figMain,'activeZ');

    cla(ax);

    FRFdata = proc.steppedsine.FRF{activeZ};
    fTest   = FRFdata.f;

    for a = 1:Na
        plot(ax, fTest, FRFdata.HMagS(a,:), 'LineWidth',1);
        hold(ax,'on');
    end

    title(ax,sprintf('FRF Magnitude — T = %d °C', temps(activeZ)));
end

% ---- drawTemperature.m ----
function drawTemperature(fig)
%DRAWTEMPERATURE Redraw backbone diagnostics for the active temperature.

    ax   = getappdata(fig,'ax');
    proc = getappdata(fig,'proc');
    Na   = getappdata(fig,'Na');
    temps = getappdata(fig,'temps');
    Z    = getappdata(fig,'activeZ');

    % --- SAVE CURRENT AXIS STATE ---
xL = ax.XLim;
yL = ax.YLim;
yScale = ax.YScale;

cla(ax);   % clear old plots

% --- RESTORE AXIS STATE ---
ax.XLim = xL;
ax.YLim = yL;
ax.YScale = yScale;
ax.XLimMode = 'manual';
ax.YLimMode = 'manual';

    hold(ax,'on'); ylim([0 100]);

    % ----- Amplitude curves (black dotted) -----
    amps = proc.steppedsine.amp{Z};

for a = 1:numel(amps)
    if isempty(amps(a).Fit.fDense), continue; end
    plot(ax, ...
        amps(a).Fit.fDense, ...
        amps(a).Fit.yFit, ...
        'k:', 'LineWidth', 0.8);
    ylim([0 100]);
end


    % ----- Backbone (solid blue) -----
    BB = proc.steppedsine.backbone{Z};
    if ~isempty(BB)
        plot(ax, BB.f, BB.y, 'b-', 'LineWidth', 3);
        ylim([0 100]);
    end

    title(ax, sprintf('Temperature = %d ^\\circC   (← / → to cycle)', temps(Z)));

    drawnow; ylim([0 100]);
end

% ---- findHalfPowerIntersection.m ----
function fI = findHalfPowerIntersection(fDense, yFit, idxPeak, direction, yHalf)
%FINDHALFPOWERINTERSECTION Interpolate a half-power crossing on a dense response curve.
    % Linear interpolation to find the exact frequency where amplitude crosses yHalf
    fI = NaN;
    n = numel(fDense);
    if strcmp(direction, 'left')
        if idxPeak <= 1, return; end
        % Look backwards from the peak for the first crossing
        idxL = find(yFit(1:idxPeak) <= yHalf, 1, 'last');
        if isempty(idxL) || idxL == idxPeak, return; end

        yseg = yFit(idxL:idxL+1);
        fseg = fDense(idxL:idxL+1);
        if diff(yseg) == 0, fI = mean(fseg); else, fI = interp1(yseg, fseg, yHalf); end
    else % 'right'
        if idxPeak >= n, return; end
        % Look forwards from the peak for the first crossing
        idxRrel = find(yFit(idxPeak:end) <= yHalf, 1, 'first');
        if isempty(idxRrel), return; end
        idxR = idxPeak - 1 + idxRrel;
        if idxR <= 1 || idxR > n || idxR-1 < 1, return; end

        yseg = yFit(idxR-1:idxR);
        fseg = fDense(idxR-1:idxR);
        if diff(yseg) == 0, fI = mean(fseg); else, fI = interp1(yseg, fseg, yHalf); end
    end
end

% ---- localPeakEstimate.m ----
function p = localPeakEstimate(vec)
%LOCALPEAKESTIMATE Return a robust local response-peak estimate for one dwell.
    % Find local maxima where a point is greater than its immediate neighbors
    p = NaN;
    if isempty(vec), return; end
    if numel(vec) < 3
        p = max(vec);
        return;
    end
    idx = find( vec(2:end-1) > vec(1:end-2) & vec(2:end-1) > vec(3:end) );
    if ~isempty(idx)
        % Return the average of the found local maxima (for robust peak estimation)
        p = mean(vec(idx+1));
    else
        p = max(vec);
    end
end

% ---- temperatureCycleCallback.m ----
function temperatureCycleCallback(fig,evt)
%TEMPERATURECYCLECALLBACK Keyboard callback for temperature-index browsing.

    activeZ = getappdata(fig,'activeZ');
    temps   = getappdata(fig,'temps');
    nT      = numel(temps);

    switch evt.Key
        case 'rightarrow'
            activeZ = min(activeZ+1,nT);
        case 'leftarrow'
            activeZ = max(activeZ-1,1);
        otherwise
            return
    end

    setappdata(fig,'activeZ',activeZ);
    drawTemperature(fig);

    figFRF = getappdata(fig,'figFRF');

if isvalid(figFRF)
    drawFRF(figFRF);
end

end

% ---- extracted from ImpactHammer2.m ----
function f_cross = findHalfPowerImpact(fD, yD, iPk, side, yH)
%FINDHALFPOWERIMPACT Interpolate the requested impact-FRF half-power crossing.
    f_cross = NaN;
    if strcmp(side,'left')
        idx = find(yD(1:iPk) <= yH, 1, 'last');
        if ~isempty(idx) && idx < iPk
            f_cross = interp1(yD(idx:idx+1), fD(idx:idx+1), yH);
        end
    else
        seg = yD(iPk:end);
        idx = find(seg <= yH, 1, 'first');
        if ~isempty(idx) && idx > 1
            i = iPk+idx-1;
            f_cross = interp1(yD(i-1:i), fD(i-1:i), yH);
        end
    end
end

% ---- extracted from SteppedSine2.m ----
function p = wpolyfit(x, y, w, deg)
%WPOLYFIT Numerically stable weighted polynomial fit with original-x coefficients.
    % Stable weighted polynomial fit, with original-x coefficients returned.
    x=x(:); y=y(:); w=w(:);
    good=isfinite(x) & isfinite(y) & isfinite(w) & w>0;
    x=x(good); y=y(good); w=w(good);
    if isempty(x)
        error('wpolyfit:NoData','No positive-weight finite observations.');
    end
    requestedDegree=deg;
    center=mean(x);
    scale=max(abs(x-center));
    if scale==0, scale=1; end
    z=(x-center)/scale;
    sw=sqrt(w/max(w));
    deg=min(deg,numel(unique(x))-1);
    while deg>=0
        V=zeros(numel(z),deg+1);
        for k=0:deg
            V(:,deg+1-k)=z.^k;
        end
        weightedV=V.*sw;
        if rank(weightedV)==deg+1, break; end
        deg=deg-1;
    end
    if deg<requestedDegree
        fprintf('  wpolyfit: degree %d reduced to %d (limited independent observations).\n',requestedDegree,deg);
    end
    coeffZ=weightedV\(sw.*y);
    p=coeffZ(1);
    for k=2:numel(coeffZ)
        p=conv(p,[1/scale,-center/scale]);
        p(end)=p(end)+coeffZ(k);
    end
    p=[zeros(1,requestedDegree-deg),p(:).'];
end

% ---- extracted from TMCoupling.m ----
function val = getVar(m, varname)
%GETVAR Return a named fitted-model coefficient, or zero when the term is absent.
    idx = strcmp(m.CoefficientNames, varname);
    if any(idx); val = m.Coefficients.Estimate(idx); else; val = 0; end
end

% ---- extracted from fitDuffingAlpha ----
function [T_vec, alpha_vec, f0sq_vec, R2_vec] = fitDuffingAlpha(obs)
%FITDUFFINGALPHA Fit f_res^2 = f0^2 + alpha*A^2 independently at each temperature.
    temps = unique(obs.Temperature_degC);
    T_vec    = NaN(size(temps));
    alpha_vec = NaN(size(temps));
    f0sq_vec  = NaN(size(temps));
    R2_vec    = NaN(size(temps));

    for ti = 1:numel(temps)
        T  = temps(ti);
        rows = obs(obs.Temperature_degC == T, :);

        % Need resonance frequency and displacement amplitude
        if ~any(strcmp(rows.Properties.VariableNames,'ResonanceFreq_Hz')) || ...
           ~any(strcmp(rows.Properties.VariableNames,'A_metres'))
            warning('PlotDuffingAlpha: missing columns in table.'); continue;
        end

        f_res = rows.ResonanceFreq_Hz;
        A_m   = rows.A_metres * 1e6;   % convert to µm

        valid = ~isnan(f_res) & ~isnan(A_m) & A_m > 0;
        if sum(valid) < 4; continue; end

        f2 = f_res(valid).^2;
        A2 = A_m(valid).^2;

        % Least-squares:  f^2 = f0^2 + alpha * A^2
        X    = [ones(sum(valid),1), A2];
        beta = X \ f2;

        f0sq  = beta(1);
        alpha = beta(2);

        f2_pred = X * beta;
        ss_res  = sum((f2 - f2_pred).^2);
        ss_tot  = sum((f2 - mean(f2)).^2);
        R2 = 0;
        if ss_tot > 0; R2 = 1 - ss_res/ss_tot; end

        T_vec(ti)     = T;
        alpha_vec(ti) = alpha;
        f0sq_vec(ti)  = f0sq;
        R2_vec(ti)    = R2;
    end
end

% ---- shared damping regression utility ----
function [log_c_ref_out, d1_out] = fitLogDamping(dT_in, c_in, d1_bound)
%FITLOGDAMPING Fit the bounded log-linear equivalent-damping temperature model.
    % Log-scale regression for c_ref and d1.
    % Clamps d1 to [-d1_bound, +d1_bound] to prevent unphysical extrapolation.
    pos = c_in > 0;
    if sum(pos) < 3
        log_c_ref_out = log(max(median(c_in),eps));
        d1_out        = 0;
        return;
    end
    Phi_log = [ones(sum(pos),1), dT_in(pos)];
    b_log   = Phi_log \ log(c_in(pos));
    log_c_ref_out = b_log(1);
    d1_out        = max(-d1_bound, min(d1_bound, b_log(2)));
end
