# vasp_emu: An emulator for machine learning potentials in VASP.

## How to install:

### Add these variables to your ~/.bashrc (or wherever you store environment variables):
```
export PATH=$PATH:/path/to/vasp_emu/bin/
export PYTHONPATH=$PYTHONPATH:/path/to/vasp_emu/
export PYTHONPATH=$PYTHONPATH:/home/ma62425/vasp_emu/tsase/
export ASE_VASP_COMMAND="mpirun path/to/vasp/executable"
```
**Acessing UMA models:**
1. Set up a HuggingFace account here: https://huggingface.co
2. Click on your account profile and select "Access Tokens"
3. Generate an access token.
4. In your terminal, type "huggingface-cli login" and enter the access token.

Now you have access to the UMA models. These instructions follow those given by the people at fairchem. Their set of instructions are here: https://fair-chem.github.io/core/install.html .

**There is many different functionalities offered in vasp_emu and an example for each:**
- ML Relaxation:
  `examples\O2_relaxation`

- ML Relaxation followed by DFT relaxation:
  `examples\ml_w_vasp`

- ML NEB:
  `examples\MeOH_rxn_neb`

- ML Molecular Dynamics:
  `examples\Cu_MD`

- ML Molecular Dynamics with DFT checks:
  coming soon
