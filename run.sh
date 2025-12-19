#!/bin/bash

### Runner script for TVGuide-Aniso
### Script will check for a local tvguide environment and the mamba/conda commands
### If not found, will initiate build script
### If found, will start tvguide-aniso with custom gunicorn settings
### Gunicorn settings can be modified in order to fit system settings as necessary


if [[ `basename "$PWD"` != "tvguide-aniso" ]]; then
    git clone https://github.com/widi9545/tvguide-aniso
    cd tvguide-aniso
fi 

export MINIFORGE_LOCATION=$(dirname ${BASH_SOURCE[0]})


bash $MINIFORGE_LOCATION/build.sh
if [[ -z `command -v conda` || -z `command -v mamba`  ]]; then
    MINIFORGE_PREFIX="$MINIFORGE_LOCATION/.miniforge3"
    if [ -d $MINIFORGE_PREFIX ]; then
        source "$MINIFORGE_PREFIX/etc/profile.d/conda.sh"
    else
        echo "Mamba/Conda commands not found - exiting"
        exit 1
    fi
fi

echo "Starting TVGuide-Aniso "
eval "$(conda shell.bash hook)"
conda activate $MINIFORGE_LOCATION/tvguide
pkill -9 gunicorn
gunicorn -w 3 -t 6 -b 128.138.136.178:8000 tvguide:server
