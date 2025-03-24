CoalPT
======

[![Build Status](https://github.com/nicfel/CoalPT/workflows/Unit%2Fintegration%20tests/badge.svg)](https://github.com/nicfel/CoalPT/actions?query=workflow%3A%22Unit%2Fintegration+tests%22)


BEAST 2 package for inference under the coalescent with plasmid transfer **CoalPt**. CoalPt models chromosomal and plasmid DNA co-evolution using a joint coalescent and plasmid transfer process in a Bayesian phylogenetic network approach. This approach reconstructs differences in the evolutionary history of plasmids and chromosomes to reconstruct instances where plasmids likely move between bacterial lineages while accounting for uncertainty in the data.

This repository contains the source code for CoalPt. The package will become part of the BEAST2 package manager upon publication. Prior to that, the easies way to install CoalPt is to use the latest release. To see where packages are installed on your computer, you can open the package manager in BEAUti and click the **?** symbol (on mac it's typically /Users/username/Library/Applications Support/BEAST/2.7). If you unpack CoalPt into that repository, it should show up in the package manager and be usable.



Building CoalRe
---------------

In order to build CoalPt from the source, you will need the following:

1. [OpenJDK](https://adoptopenjdk.net) v8 or later,
2. The Apache Ant build tool.

Once these are installed, open a shell in the root directory of this repository
and use

    $ ant

to build the package.

You might need to include

`JAVA_FX_HOME=/Library/Java/JavaVirtualMachines/zulu-17.jdk/Contents/Home/`

License
-------

CoalRe is free (as in freedom) software and is distributed under the terms of
version 3 of the GNU General Public License.  A copy of this license is found
in the file named `COPYING`.
