# EVs_transport
In this repository, you can find the source code used to track extracellular vesicles (EVs) motion from microscopy videos and simulate in-silico the prion protein-mediated transport of EVs on the neuron surface.
The model and the results of the numerical simulations are contained in this paper currently under review [[1]](#1).

## Repository structure
The repository is structured as follows:

* `tracking`: contains : 
	* `Detect_positions_MVesicle.m`: the code for extracting the EVs position from the microscopy movies in [[2]](#2);
	* `coord_processing.m`: the code for processing the txt files obtained from `Detect_positions_MVesicle.m`;

* `simulations`: contains :
	* `EVs_model.m`: the code for the numerical simulations of the mathematical model proposed in [[1]](#1).

## Citing

If you find this code useful for your work, please cite [[1]](#1).

## Licence

The source code contained in this repository is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 2.1 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

## References

<a id="1">[1]</a>
G. Pozzi, G. Mazzilli, G. D'Arrigo, C. Verderio, G. Legname, S. Turzi, P. Ciarletta (2025).
Modeling the prion protein-mediated transport of extracellular vesicles on the neuron surface.
arXiv preprint arXiv:2502.03610, under review.

<a id="1">[2]</a>
G. Pozzi, G. Mazzilli, G. D'Arrigo, C. Verderio, G. Legname, S. Turzi, P. Ciarletta (2025).
Video Modeling the prion protein-mediated transport of extracellular vesicles on the neuron surface.
10.5281/zenodo.15031251
