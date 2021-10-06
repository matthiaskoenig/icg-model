# ICG model
This repository provides the indocyanine green physiological based pharmacokinetics model (PBPK) described in

> Prediction of survival after partial hepatectomy using a physiologically based pharmacokinetic model of indocyanine green liver function tests
*Adrian Köller, Jan Grzegorzewski, Hans-Michael Tautenhahn, Matthias König*
bioRxiv 2021.06.15.448411; doi: [https://doi.org/10.1101/2021.06.15.448411](https://doi.org/10.1101/2021.06.15.448411)


The model is distributed as [SBML](http://sbml.org) and available from [`icg_body_flat.xml`](./models/icg_body_flat.xml)

![Visualization of the model using cy3sbml](./visualization/icg_body_flat.png)

The SBML4humans model report is available via
[https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/icg-model/main/models/icg_body.xml](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/icg-model/main/models/icg_body_flat.xml)

### Comp submodels
The liver submodel is available from [`icg_liver.xml`](./models/icg_liver.xml) with corresponding SBML4humans report at
[https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/icg-model/main/models/icg_body.xml](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/icg-model/main/models/icg_liver.xml)

![Visualization of the liver submodel using cy3sbml](./visualization/icg_liver.png)

The whole-body submodel is available from [`icg_body.xml`](./models/icg_body.xml) with corresponding SBML4humans report at
[https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/icg-model/main/models/icg_body.xml](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/icg-model/main/models/icg_body.xml)


## How to cite

 Prediction of survival after partial hepatectomy using a physiologically based pharmacokinetic model of indocyanine green liver function tests
*Adrian Köller, Jan Grzegorzewski, Hans-Michael Tautenhahn, Matthias König*
bioRxiv 2021.06.15.448411; doi: [https://doi.org/10.1101/2021.06.15.448411](https://doi.org/10.1101/2021.06.15.448411)


## License

* Source Code: [LGPLv3](http://opensource.org/licenses/LGPL-3.0)
* Documentation: [CC BY-SA 4.0](http://creativecommons.org/licenses/by-sa/4.0/)
* Models: [CC BY-SA 4.0](http://creativecommons.org/licenses/by-sa/4.0/)

The icg-model source is released under both the GPL and LGPL licenses version 2 or
later. You may choose which license you choose to use the software under.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License or the GNU Lesser General Public
License as published by the Free Software Foundation, either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

Funding
=======
Adrian Köller and Matthias König are supported by the Federal Ministry of Education and Research (BMBF, Germany)
within the research network Systems Medicine of the Liver (**LiSyM**, grant number 031L0054). Matthias König
is supported by the German Research Foundation (DFG) within the Research Unit Programme FOR 5151
[QuaLiPerF](https://qualiperf.de) (Quantifying Liver Perfusion-Function Relationship in Complex Resecti)
A Systems Medicine Approach)" by grant number 436883643.

© 2020-2021 [Matthias König](https://livermetabolism.com)
