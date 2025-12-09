# vasp_emu: An emulator for machine learning potentials in VASP.

## How to install:

### Add the /bin directory to your PATH variable and path/to/vasp_emu/ to the $PYTHONPATH variable:
```
export PATH=$PATH:/path/to/vasp_emu/bin/
export PYTHONPATH=$PYTHONPATH:/path/to/vasp_emu/
```
Acessing UMA models:
Step 1: Set up a HuggingFace account here: https://huggingface.co
Step 2: Click on your account profile and select "Access Tokens"
Step 3: Generate an access token.
Step 4: In your terminal, type "huggingface-cli login" and enter the access token.

Now you have access to the UMA models. These instructions follow those given by the people at fairchem. Here is are those set of instructions: https://fair-chem.github.io/core/install.html .

Head over to \examples\O2_relaxation for a quick start with vasp_emu.
