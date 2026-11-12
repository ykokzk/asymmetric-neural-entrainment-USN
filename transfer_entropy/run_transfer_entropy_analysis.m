% run_transfer_entropy_analysis.m
%
% Transfer entropy (TE) analysis used in Okazaki et al., "Asymmetric neural
% entrainment at resonance frequencies underlies unilateral spatial neglect".
%
% This script computes the transfer entropy between the amplitude envelopes of
% band-pass filtered EEG signals, for every channel pair and every delay.
% It covers the TE computation of the original analysis only. Selection of the
% 9-Hz flicker condition, region-of-interest selection, taking the maximum over
% delays, group statistics and figure generation were done in separate
% downstream code and are NOT part of this script.
%
% The preprocessed EEG data used in the study are not distributed. To run this
% script, load your own preprocessed data into the workspace in the format
% described below, then run this script.
%
% ---------------------------------------------------------------------------
% REQUIRED INPUT (must already exist in the workspace)
%
%   data.segdata.<condition> : [channels x time x trials]
%       One field per experimental condition. Every field is processed, and the
%       fields form the 4th dimension of the output in the order returned by
%       fieldnames(data.segdata). The first 19 channels are analysed, in the
%       channel order of the recording montage (the original study used a
%       19-channel 10-20 montage). Sampling frequency: 500 Hz.
%
%   data.timeinfo : [1 x time]
%       Time axis in seconds, with the same length as dimension 2 of segdata
%       and 0 s at stimulus onset. Used to cut the analysed segment
%       (-10 to 20 s) out of each trial.
%
%   The band-pass filter and the Hilbert envelope are applied by this script,
%   so data.segdata must contain broadband (not narrow-band filtered) EEG.
%
% ---------------------------------------------------------------------------
% DEPENDENCIES
%
%   MATLAB (base)
%   Signal Processing Toolbox
%       hilbert : amplitude envelope
%       fir1    : used internally by the FieldTrip FIR band-pass filter
%   FieldTrip (not distributed here; add it to your MATLAB path)
%       ft_preproc_bandpassfilter : 8.5-9.5 Hz FIR band-pass filter
%       nearest                   : converts the segment limits (s) to samples
%   TE.m (distributed alongside this script)
%       Transfer entropy class provided by H. Yokoyama.
%
% ---------------------------------------------------------------------------
% OUTPUT (workspace variable output_data_te)
%
%   output_data_te.TEy2x : [epoch x channel pair x delay x condition]
%       Transfer entropy from Y to X. For channel pair k, X is channel
%       chcomb(k,1) and Y is channel chcomb(k,2).
%   output_data_te.TEx2y : [epoch x channel pair x delay x condition]
%       Transfer entropy from X to Y, with the same channel-pair convention.
%
%       Both are averaged over trials, as in the original analysis; single
%       trials are not kept. Dimension 1 is the epoch (with the settings below:
%       1 = -10 to 0 s, 2 = 0 to 10 s, 3 = 10 to 20 s), NOT the trial.
%
%   output_data_te.dimord         : 'epc_chcmb_testep_cond' (as in the original analysis)
%   output_data_te.chcomb         : [channel pairs x 2], nchoosek(1:19,2)
%   output_data_te.delay_sample   : [1 x delays] delay of each delay index (samples)
%   output_data_te.delay_sec      : [1 x delays] delay of each delay index (s)
%   output_data_te.epoch          : [epochs x 2] first and last sample of each epoch
%   output_data_te.condition_list : condition names, in the order of dimension 4
%   output_data_te.timeinfo       : time axis (s) of the analysed segment
%
% ---------------------------------------------------------------------------

%% Analysis parameters (values used for the published analysis)
sf = 500;               % sampling frequency (Hz)

tewin  = 5;             % TE window (s)
testep = 0.1;           % TE delay step (s)

savesegment = [-10 20]; % analysed segment relative to stimulus onset (s)

epochsize = 10;         % epoch length (s)
epochstep = epochsize;  % epoch step (s): non-overlapping epochs

doFilt     = 1;         % apply the FIR band-pass filter
doEnvelope = 1;         % take the Hilbert amplitude envelope
Fbp        = [8.5 9.5]; % band-pass filter band (Hz)

%% Check that the input data are in the workspace
if ~exist('data','var') || ~isstruct(data) || ~isfield(data,'segdata') || ~isfield(data,'timeinfo')
    error(['Load your preprocessed EEG data into the workspace before running this script. ', ...
           'Required: data.segdata.<condition> [channels x time x trials] and ', ...
           'data.timeinfo [1 x time] in seconds.']);
end

%% Convert seconds to samples
start_point_save = nearest(data.timeinfo, savesegment(1));
end_point_save   = nearest(data.timeinfo, savesegment(2));
testep      = round(testep*sf);
tewinsample = round(tewin*sf);

%% Main loop
output_data_te = [];
output_data_te.timeinfo = data.timeinfo(start_point_save:end_point_save);

condition_list = fieldnames(data.segdata);
output_data_te.condition_list = condition_list;
chcomb = nchoosek(1:19,2);
output_data_te.chcomb = chcomb;

for c = 1:length(condition_list)  % loop for conditions
    fprintf('Processing condition %d/%d: %s\n', c, length(condition_list), condition_list{c});

    input_data = getfield(data.segdata, condition_list{c});
    input_data = input_data(1:19,:,:);

    tmp2TEy2x = []; tmp2TEx2y = [];
    epcstepsample = epochstep*sf;
    epcsizesample = epochsize*sf;
    for tri = 1:size(input_data,3) % loop for trial
        tmpinput = input_data(1:19,start_point_save:end_point_save,tri);

        % bandpass filter
        if doFilt
            [tmpinput] = ft_preproc_bandpassfilter(tmpinput, sf, Fbp, [], 'fir'); % fir=fir1
            if doEnvelope
                y = hilbert(tmpinput'); % transpose ch x time -> time x ch, otherwise hilbert is applied along the wrong dimension
                tmpinput = abs(y');
            end
        end

        nepoch  = round(size(tmpinput,2)/epcsizesample); % with epochsize=10, the -10 to 20 s data are split into 3 epochs
        epc     = 0:nepoch-1;
        nchcomb = length(chcomb);
        ndelay  = length(2:testep:epcsizesample-tewinsample);
        tmpTEy2x = zeros(nepoch,nchcomb,ndelay);
        tmpTEx2y = zeros(nepoch,nchcomb,ndelay);

        for epci = 1:length(epc) % loop for pre, during, post epochs

            epochtmp   = [epc(epci)*epcstepsample+1 (epc(epci)+1)*epcstepsample];
            tmpwindata = tmpinput(:,epochtmp(1):epochtmp(2)); % pre, during, post data

            testy2x = zeros(nchcomb,ndelay);
            testx2y = zeros(nchcomb,ndelay);
            for tmpch = 1:nchcomb                % channel combination
                X = tmpwindata(chcomb(tmpch,1),1:tewinsample);
                Y = tmpwindata(chcomb(tmpch,2),1:tewinsample);
                teti = 0;
                for tet = 2:testep:epcsizesample-tewinsample   % compute for each TE delay
                    teti = teti + 1;
                    Xdelay = tmpwindata(chcomb(tmpch,1),tet:tet+tewinsample-1); % series shifted by tet relative to X
                    Ydelay = tmpwindata(chcomb(tmpch,2),tet:tet+tewinsample-1);
                    te = TE(X,Y,Xdelay, Ydelay);
                    te.get_PDF()
                    [testy2x(tmpch,teti), testx2y(tmpch,teti)] = te.get_TE();
                end
            end

            tmpTEy2x(epci,:,:) = testy2x;
            tmpTEx2y(epci,:,:) = testx2y;
        end

        tmp2TEy2x(:,:,:,tri) = tmpTEy2x;
        tmp2TEx2y(:,:,:,tri) = tmpTEx2y;
    end

    output_data_te.TEy2x(:,:,:,c) = squeeze(mean(tmp2TEy2x,4)); % average trial
    output_data_te.TEx2y(:,:,:,c) = squeeze(mean(tmp2TEx2y,4)); % average trial
end

output_data_te.epoch  = [1 epcsizesample]+epcstepsample*[0 1 2]';
output_data_te.dimord = 'epc_chcmb_testep_cond';

% Delay axis of dimension 3 of TEy2x / TEx2y. The delayed series starts at
% sample tet of the epoch while the reference series starts at sample 1, so the
% delay is tet-1 samples.
delay_startsample = 2:testep:epcsizesample-tewinsample;
output_data_te.delay_sample = delay_startsample - 1;
output_data_te.delay_sec    = (delay_startsample - 1)/sf;
