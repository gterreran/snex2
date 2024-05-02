#!/bin/bash
# file: create_snex2_env.sh

echo "Building the SNEx2 data folder tree."
echo

while true ; do
    read -r -p "Enter where you want the SNEx2 data (including the database) to be stored [~]: " filepath
    filepath=${filepath:-~}
    filepath="`eval echo ${filepath//>}`"
    if [ -d $filepath ] ; then
        break
    fi
    echo "$filepath not found. Please insert a correct path"
done


folders=("$filepath/supernova" "$filepath/supernova/snex_images" "$filepath/supernova/snex_images/thumbs" "$filepath/supernova/snex2" "$filepath/supernova/snex2/db_data" "$filepath/supernova/data" "$filepath/supernova/data/lsc" "$filepath/supernova/data/floyds" "$filepath/supernova/data/gw" "$filepath/supernova/data/2m0a" "$filepath/supernova/data/fts" "$filepath/supernova/data/0m4")
ask=true

for d in "${folders[@]}" ; do
    if [ ! -d $d ] ; then
        mkdir $d
    else
        #If folder exists ask if the user wants to erase it
        while true; do
            #If $ask is false, it won't ask again
            if [ "$ask" = true ]; then
                read -p "Folder $d already exists. Do you want to erase it (recomanded for a clean installation) ['y,n,Y,N] " yn
            fi
            case $yn in
                [Yy]* ) rm -r $d; mkdir $d; break;; #erasing the folder and recreating it
                [Nn]* ) echo "Existing folder $d kept as it is."; break;; #Not erising it
                * ) echo "Please answer yes or no.";;
            esac
        done
        #If the user replied Y or N, set $ask to false so it won't ask again
        if [[ $yn =~ ^[YN]$ ]]; then
                ask=false
        fi
    

    fi

done

printf "\n\nThe SNEx2 data folders are now organized as the following tree\n\n"
printf ${filepath}"/supernova\n"
printf "|-- snex_images\n"
printf "|   \`-- thumbs\n"
printf "|-- snex2\n"
printf "|   \`-- db_data\n"
printf "|-- data\n"
printf "    |-- lsc\n"
printf "    |-- floyds\n"
printf "    |-- gw\n"
printf "    |-- 2m0a\n"
printf "    |-- fts\n"
printf "    \`-- 0m4\n"

echo

echo "Creating the yaml file."
echo

#Considering that we create the yaml file from scratch, we can even consider
#to prompt the user to insert their own user and password for the database

# cat > snex2_conda_environment.yml << EOF
# # run: conda env create --file snex2_conda_environment.yml"
# name: snex2
# dependencies:
#  - python=3.10
#  - mysql
#  - postgresql
#  - conda-forge::ligo.skymap
#  - pip
#  - pip:
#    - -r requirements.txt
# variables:
#  SNEX2_DB_BACKEND: postgres
#  SNEX2_DB_USER: snex2dbuser
#  SNEX2_DB_PASSWORD: snex2dbpassword
#  SNEX2_DB_HOST: 127.0.0.1
#  SNEX2_DB_PORT: 5435
#  SNEX2_DB_DATA_PATH: $filepath/supernova/snex2/db_data
# EOF

conda create -n snex2 -c conda-forge python=3.10 mysql postgresql ligo.skymap
conda env config vars set SNEX2_DB_BACKEND=postgres SNEX2_DB_USER=snex2dbuser SNEX2_DB_PASSWORD=snex2dbpassword SNEX2_DB_HOST=127.0.0.1 SNEX2_DB_PORT=5435 SNEX2_DB_DATA_PATH=$filepath/supernova/snex2/db_data --name snex2
eval "$(conda shell.bash hook)" #necessary to activate the environment
conda activate snex2
pip install -r requirements.txt
conda deactivate

printf "The SNEx2 environment is now set. To activate it, type\n\n  conda activate snex2\n\n
echo "DONE!"



