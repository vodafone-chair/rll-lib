# Readme
MATLAB library for runlength-limited (RLL) coding. This library contains encoders and decoders, to utilize the RLL codes derived in [1], [2]. The functions can, e.g., be used to implement zero-crossing modulation (ZXM).

## Getting Started
The file `minimumWorkingExample.m` provides a quick start.

## List of Main Functions
- `rllEncoderFct.m` performs RLL encoding using the FSM RLL codes derived in [1], [2]
- `rllBcjrFct.m` performs RLL soft-input soft-output decoding by implementing the Bahl-Cocke-Jelinek-Raviv (BCJR) algorithm as proposed in [1], [2]
- `rllViterbiFct.m` performs minimum Hamming distance decoding of the RLL sequences by implementing the Viterbi algorithm as proposed in [Algorithm 1, 4]
- `getFsmRllAutocorrFct.m` computes the autocorrelation function of the RLL sequences generated by the FSM encoders as described in [1], [3]

## References
[1] P. Neuhaus, M. Dörpinghaus, and G. Fettweis, “Zero-Crossing Modulation for Wideband Systems Employing 1-Bit Quantization and Temporal Oversampling: Transceiver Design and Performance Evaluation,” in IEEE Open Journal of the Communications Society, vol. 2, pp. 1915-1934, 2021, doi: 10.1109/OJCOMS.2021.3094927.\
[2] P. Neuhaus, M. Dörpinghaus, H. Halbauer, S. Wesemann, M. Schlüter, F. Gast, and G. Fettweis, “Sub-THz  wideband  system  employing  1-bit  quantization  and  temporal  oversampling,” in Proc. IEEE Int. Conf. Commun. (ICC), Dublin, Ireland, pp. 1–7, Jun. 2020, doi: 10.1109/ICC40277.2020.9148753.\
[3] P. Neuhaus, M. Dörpinghaus, H. Halbauer, V. Braun, and G. Fettweis, “On the spectral efficiency of oversampled 1-bit quantized systems for wideband LOS channels,” in Proc. IEEE Int. Symp. on Personal, Indoor and Mobile Radio Commun. (PIMRC), London, U.K., pp. 1-6, Aug. 2020, doi: 10.1109/PIMRC48278.2020.9217277.\
[4] P. Neuhaus, D. M. V. Melo, L. T. N. Landau, R. C. de Lamare, and G. Fettweis. “Zero-Crossing Modulations for a Multi-User MIMO Downlink with 1-Bit Temporal Oversampling ADCs”. In Proc. European Signal Proc. Conf. (EUSIPCO). Dublin, Ireland, pp. 816-820, Sep. 2021, doi: 10.23919/EUSIPCO54536.2021.9616302.
