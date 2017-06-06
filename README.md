# Silent Substitution Toolbox

#### 1. Description
MATLAB toolbox to compute estimates of human photoreceptor spectral sensitivities, compute silent substitution spectral modulations for a variety of devices (including standard monitors and devices with many narrowband primaries), and to estimate the normal variation of photoreceptor spectral sensitivities and how much contrast splatter a given modulation can be expected to produce on nominally silenced photoreceptors.

#### 2. License
This software is licensed under the license specified in `LICENSE.md` (MIT License). The code may be used freely. Using it to compute robust modulations by silencing multiple photoreceptor classes is covered by a U.S. Patent Application ([United States Patent Application 14/852001, "ROBUST TARGETING OF PHOTOSENSITIVE MOLECULES", 09/11/2015](http://www.freepatentsonline.com/y2016/0073922.html)) and is subject to licensing.

#### 3. Developers
This software was developed by: 
* [Manuel Spitschan](https://github.com/spitschan) (Lead Developer; Stanford University, 2016-now, University of Pennsylvania, 2012-2016)
* [Geoffrey K. Aguirre](https://github.com/gkaguirre) (University of Pennsylvania)
* [David H. Brainard](https://github.com/DavidBrainard) (University of Pennsylvania)

#### 4. Citing
If you use this code in support of work in a published paper, please cite us:

* Spitschan M, Aguirre GK, Brainard DH (2015) Selective Stimulation of Penumbral Cones Reveals Perception in the Shadow of Retinal Blood Vessels. _PLoS ONE 10_(4): e0124328. [doi: 10.1371/journal.pone.0124328](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0124328)

#### 5. Dependencies and pre-requisites
The *Silent Substitution Toolbox* relies on functions from [Psychtoolbox-3](https://github.com/Psychtoolbox-3/Psychtoolbox-3) and the [Brainard Lab Toolbox](https://github.com/BrainardLab/BrainardLabToolbox). These dependencies can be obtained manually and put on the path. Alternatively, they can be obtained using [ToolboxToolbox](https://github.com/ToolboxHub/ToolboxToolbox) (TbTb). Once TbTb is installed, simply run:

```
tbUse('SilentSubstitutionToolbox');
```

#### 6. Use
For an overview of the functions contained within this toolbox, please look at `Contents.m`.

#### 7. Version history

* [v1.2.1](https://github.com/spitschan/SilentSubstitutionToolbox/releases/tag/v1.2.1) (April 24, 2017)
* [v1.2](https://github.com/spitschan/SilentSubstitutionToolbox/releases/tag/v1.2) (January 15, 2017)
* [v1.1](https://github.com/spitschan/SilentSubstitutionToolbox/releases/tag/v1.1) (January 29, 2015)
* [v1.0](https://github.com/spitschan/SilentSubstitutionToolbox/releases/tag/v1.0) (December 14, 2014)

#### 7. Questions, Feedback?
Please open an issue on the GitHub page for this toolbox. For other questions, [email Manuel Spitschan](mailto:spitschan@stanford.edu).
