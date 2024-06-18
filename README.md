![workflow badge](https://github.com/Omer-Sella/802.3/actions/workflows/bchEncoder.yml/badge.svg)
![workflow badge](https://github.com/Omer-Sella/802.3/actions/workflows/clause177_5_encoder.yml/badge.svg)

## BCH and Reed Solomon (RS) decoding
The decoder for BCH and RS codes is from this [repository](https://github.com/Omer-Sella/reedSolomon "Reed Solomon and BCH decoding"). It also has an encoder (which is just polynomial division), which I used to check encoder definition in the 802.3 physical layer amendment.

## Clause 177.5 decoding
A chase-2 decoder that is riding on a hard decision Hamming decoder.


## Bucket
Clause 177.5 using Viterbi that relies on a Viterbi decoder from my [Viterbi repository](https://github.com/Omer-Sella/viterbi "Viterbi decoder implementation(s)" and other work from [Euclid](https://github.com/Omer-Sella/Euclid "Encoding and decoding data into DNA"), which is where I keep most of my DNA-data-storage work.
