##depends:conda

CONDA_PATH=("/opt/chipster/tools/miniconda3")
export PATH=${PATH}:${CONDA_PATH}/bin
source activate chipster_tools
conda install minimap2=2.9 -y 