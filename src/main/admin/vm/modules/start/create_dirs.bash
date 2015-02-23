##depends:none

if [ ! -d /mnt/tools ]; then
	mkdir /mnt/tools
fi

# Symlink to tools
ln -s /mnt/tools ${TOOLS_PATH}