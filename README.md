# BScMultiHiggsPySecDec
Generation of the amplitude form factors of single and double Higgs production in gluon fusion with Mathematica, and subsequent numerical evaluation of these form factors using pySecDec.

Required software:
- Wolfram Mathematica
  * [alibrary](https://github.com/magv/alibrary)
  * [hepware](https://github.com/magv/hepware) to install the following, requisite packages:
    * [Qgraf](http://cfif.ist.utl.pt/~paulo/qgraf.html)
    * [Feynson](https://github.com/magv/feynson)
    * [FORM](https://github.com/vermaseren/form)
    * [Kira](https://gitlab.com/kira-pyred/kira)
- Python
  * [pySecDec](https://secdec.readthedocs.io/en/stable/)
  * [Jupyter](https://jupyter.org/)

In order to use the Mathematica files, you need to replace the placeholder strings or file paths `"PathToAlibrary"` with your respective path to the directory that contains your copy of _alibrary_, and `"PathToOutputDirectory[i]"` with the path to the directories you want the output `compile.py` files to be written to. For single Higgs, the placeholder string does not contain the `"[i]"` at the end, and in the case of double Higgs, this represents both output paths.
You also need to add the model file provided here to your _alibrary_ directory, or a corresponding model file of your own choice.
