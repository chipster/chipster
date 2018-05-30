##depends:none

#Miniconda3 installation:
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
CONDA_PATH=("/opt/chipster/tools/miniconda3")
bash Miniconda3-latest-Linux-x86_64.sh -b -p $CONDA_PATH 
export PATH=${PATH}:${CONDA_PATH}/bin
conda config --add channels conda-forge
conda config --add channels bioconda
conda create -n chipster_tools -y

cat << EOF > conda_execute
#!/bin/bash

conda_path=("$CONDA_PATH")

c_env=\$(echo \$1 | awk -F "/" '{print \$1}')
c_tool=\$(echo \$1 | awk -F "/" '{print \$2}')
shift
export PATH=\${PATH}:\${conda_path}/bin
source activate \$c_env
\$c_tool \$@
EOF
