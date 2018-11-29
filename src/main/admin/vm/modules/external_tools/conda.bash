##depends:none

#Miniconda3 installation:
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
CONDA_PATH=("/opt/chipster/tools/miniconda3")
bash Miniconda3-latest-Linux-x86_64.sh -b -p $CONDA_PATH 
export PATH=${PATH}:${CONDA_PATH}/bin
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

mv $HOME/.condarc $CONDA_PATH
conda create -n chipster_tools -y

cat << EOF > $CONDA_PATH/conda_execute
#!/bin/bash

conda_path=("$CONDA_PATH")

c_env=\$(echo \$1 | awk -F "/" '{print \$1}')
c_tool=\$(echo \$1 | awk -F "/" '{print \$2}')
shift
export PATH=\${PATH}:\${conda_path}/bin
source activate \$c_env
\$c_tool \$@
EOF

chmod u+x $CONDA_PATH/conda_execute

# Activate conda
#CONDA_PATH=("/opt/chipster/tools/miniconda3")
#export PATH=${PATH}:${CONDA_PATH}/bin
source activate chipster_tools

# minimap
conda install minimap2=2.9 -y

# weasyprint
conda install -c jlmenut weasyprint -y
pip install -U html5lib=="0.9999999"