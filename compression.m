function [sparse_rep, bpp] = sparsify_W1_mse_v2(Data,cr,wname)
% This function sparsifies the rawData looking at each sensor over time
% individually (1D). It sparsifies by taking the wavelet transform and then 
% quantizing.
% now stopping based on the compression ratio

% calculate MSE and return energy ratio..

% INPUTS:
    % rawData = data to be sparsified; - [(rows x col) x time]
    % nnz_indx = indexes with values
    % nmse = acceptable error
    % wname = wavelet used; example = 'sym3'
% OUTPUTS:
    % sparse_rep = sparse representation after quantization, de-meaned
    % sparse_recon = reconstruction using sparse_rep
    % book_keeping = book keeping vector for wavelet x-form
    % quant = quantiziation value --> divide by this number and round
    %                                                       towards 0
    % nmses = evolution of nmse during search
    % quants = evolution of quants during binary search
    % sparsity = sparsity (nnz / original size)
    % q_max = max value after quantization
    % bpp = average # of bits per pixel
    % energy_ratio = energy ratio between sparse_rep and original wavelet x-form
    % means = average value for each sensor

%% Run this cell and uncomme
% nt below for example
% initialize
raw = Data;
[l,t] = size(Data);

% if nmse == 0
%     tolerance = 1;
% else
%     tolerance = strsplit(num2str(nmse),'.');
%     tolerance = numel(tolerance{2});
%     tolerance = 0.5*10^-(tolerance+1);
% end %tolerance is one extra decimal point past mse

tolerance = 0.1*cr;

% De-mean rawData --> this makes energy conservation more informative
means = zeros(l);
for row = 1:l
    means(row) = mean(Data(row,:));
    Data(row,:) = Data(row,:) - mean(Data(row,:));
end

nnz_raw = raw;

% Take transform
% rawData = [rows x col x time];
[wc,book_keeping] = wavedec(squeeze(Data(1,:)),fix(log2(t)),wname);
disp(['This is wc: ' num2str(wc)]);
% do the first one to get the size
new_sz = length(wc);
sparse_rep = zeros(l,new_sz);
for row = 1:l
    x = squeeze(Data(row,:));
    [wc,book_keeping] = wavedec(x,fix(log2(t)),wname);
    sparse_rep(row,:) = wc;
end

original_energy = norm(sparse_rep(:)); % for energy_ratio

% Start values for binary search
mx = max(sparse_rep,[],'all'); 
mn = 1e-20;
quant = mx;

% Search for best quantizer
searching = 1;


iterations = 0;

while (searching)
    % Quantize
    codeword = fix(sparse_rep./quant);
    disp(['This is a codeword: ' num2str(codeword)]);
    disp(['This is my quant: ' num2str(quant)])
    temp_sp = quant*(codeword);
    disp(['iteration number: ' num2str(iterations) ' temp_sp: ' num2str(temp_sp)]);
    % calculating the compression ratio
    % what is the maximum quantized value?
    old_quant = quant;
    q_max = round(max(temp_sp./quant,[],'all'));
    num_nnz_bits = nnz(temp_sp) * ceil((log2(q_max))+1);
    bpp = num_nnz_bits / length(temp_sp(:));

    if bpp > cr % not enough compression, quantize more
        mn = quant;
        quant = (quant + mx)/2;
    end
    if bpp < cr % too much compression, quantize less
        mx = quant;
        quant = (quant + mn)/2;
    end

    if (cr - tolerance < bpp) && (bpp < cr + tolerance) % if correct value is reached
        searching = 0;
    end



    iterations = iterations + 1;

    % if search fails, if number of iterations > 50, stop searching
    if iterations > 50 % if values stop changing
        searching = 0;
        % quant = 1e10;
        % NMSE = 1e10;
    end

end


% Estimate sparsity
sparse_rep = temp_sp;
sparsity = nnz(sparse_rep) / (length(sparse_rep(:)));
% calculate number of nnz coefficients

% what is the maximum quantized value?
q_max = round(max(sparse_rep./quant,[],'all'));
num_nnz_bits = nnz(sparse_rep) * ceil((log2(q_max))+1);
bpp = num_nnz_bits / length(sparse_rep(:));

% energy conservation
comp_power = norm(sparse_rep(:));
energy_ratio = comp_power / original_energy; %compressed power over original 
%{
disp(['sparse rep: ', num2str(sparse_rep)]);
disp(['final quant: ', num2str(quant)]);
disp(['codewords: ', num2str(sparse_rep./quant)]);
disp(['Sparsity: ' num2str(sparsity)]);
disp(['bpp: ' num2str(bpp)]);
disp(['#nnz: ' num2str(nnz(sparse_rep))]);
disp(['q_max: ' num2str(q_max)]);
disp(['Bits: ' num2str(num_nnz_bits)]);
disp(['||comp||: ' num2str(comp_power)]);
disp(['Energy ratio: ' num2str(energy_ratio)]);
%}
end