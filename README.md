# XMG - eXecutable Model Generator

XMG is a tool for the modelling of signalling pathways in Systems Biology. It is implemented in Python3 and its code underlies the PEP8 guidelines.

XMG converts [SBML](http://www.sbml.org/) models with [CellDesigner](http://www.celldesigner.org/) annotations to [SBtab](http://www.sbtab.net/) files, which directly correspond to [rxncon](http://www.rxncon.org/) input. If any of these frameworks and concepts is enigmatic for you, please follow the links to further the understanding of this tool.

While the SBtab format simply serves as intermediate format of this tool, the rxncon output file can directly be converted to an executable model. Thus, the circle from an SBML-CD map to an executable model is closed.

This repository holds a directory with development files, which were required for working on the XMG concept (executable_model_generator/concept_development).

Furthermore, the developmental status of the tool is located in directory executable_model_generator/standalone_version).
