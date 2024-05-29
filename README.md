<h1 align="center">
molvizpy
</h1>

<br>


{'MolVizPy': 'A Python package for molecular visualization, symmetry analysis. Easily visualize chemical structures in both 2D and 3D and determine point groups'}

## Usage

```python
# Install the package the the desired (example to \documents)
git clone https://github.com/alexlauris/molvizpy.git
cd \path\to\documents
cd molvizpy
cd molvizpy # Yes two times
conda env create -f molvizpy.yml
streamlit run molvizpy.py
# To stop the program, press Ctrl + C in the terminal
```

This usage example shows how to quickly leverage the package's main functionality with just one line of code (or a few lines of code). 
After importing the `main_func` (to be renamed by you), you simply pass in your `data` and get the `result` (this is just an example, your package might have other inputs and outputs). 
Short and sweet, but the real power lies in the detailed documentation.

## Installation

Create a new environment, you may also give the environment a different name! 

```
conda create -n molvizpy python=3.10 
```

```
conda activate molvizpy
```

If you need jupyter lab, install it 

```
(molvizpy) $ pip install jupyterlab
```


## Development installation

Initialize Git (only for the first time). 

Note: You should have create an empty repository on `https://github.com:alexlauris/molvizpy`.

```
git init
git add * 
git add .*
git commit -m "Initial commit" 
git branch -M main
git remote add origin git@github.com:alexlauris/molvizpy.git 
git push -u origin main
```

Then add and commit changes as usual. 

To install the package, run

```
(molvizpy) $ pip install -e ".[test,doc]"
```

### Run tests and coverage

```
(conda_env) $ pip install tox
(conda_env) $ tox
```

### License
This project is licensed under the MIT License - see the LICENSE file for details.

### Contact
For any issues, questions, or suggestions, please open an issue on GitHub or contact the maintainer at alexandre.lauris@epfl.ch


