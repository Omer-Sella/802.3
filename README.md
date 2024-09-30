![workflow badge](https://github.com/Omer-Sella/802.3/actions/workflows/bchEncoder.yml/badge.svg)
![workflow badge](https://github.com/Omer-Sella/802.3/actions/workflows/clause177_5_encoder.yml/badge.svg)

# Table of contents
1. [What this repository is about](#about)
2. [BCH and Reed Solomonlementation notes](#bchrs)
3. [Chase-2 decoding](#chase2)
4. [Getting up and running](#upandrunning)
5. [License](#license)
6. [Bucket](#bucket)

## What this repository is about<a name="about"></a>
This repository follows the development of Ethernet physical layer work groups 802.3df and 802.3dj. The idea is to provide open source reference implementation, that could also be used for validation. Finally, the intent is to be used for simulation, to help decision making at the work group.

## BCH and Reed Solomon (RS) decoding<a name="bchrs"></a>
The decoder for BCH and RS codes are from this [repository](https://github.com/Omer-Sella/reedSolomon "Reed Solomon and BCH decoding"), where you can find an encoder (which is just polynomial division), which I used to check encoder definition in the 802.3 physical layer amendment.

## Chase-2 decoding<a name="chase2"></a>
A chase-2 decoder for clause 177.5 that is riding on a hard decision Hamming decoder.

## Getting up and running <a name="upandrunning"></a>
Other than cloning this project and setting an environment that has all the requirements (listed in requirements.txt), you would need to let python know where the files are. This can be done either by setting a system environment that has the path to the project, and which the modules are set to look for (REEDSOLOMON), or, hard coding it into the modules.

## License<a name="license"></a>
Feel free to use this project under the MIT license. If you decide to reference it, please use: \
@misc{ieee8023python,\
  title        = "Reference implementation of coding methods for IEEE 802.3 .df and .dj standards",\
  author       = "{Omer S. Sella}",\
  howpublished = "\url{https://github.com/Omer-Sella/802.3}",\
  year         = 2024,\
}
## Bucket
Clause 177.5 using Viterbi that relies on a Viterbi decoder from my [Viterbi repository](https://github.com/Omer-Sella/viterbi "Viterbi decoder implementation(s)" and other work from [Euclid](https://github.com/Omer-Sella/Euclid "Encoding and decoding data into DNA"), which is where I keep most of my DNA-data-storage work.
