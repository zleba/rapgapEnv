# Installation of the RAPGAP

Commands to install Rapgap with it's dependencies (hepMC2, LHAPDF5, Pythia6).
See the content of the `Dockerfile`.

## Running it

Clone the repository to your machine
``
git clone https://github.com/zleba/rapgapEnv.git
cd rapgapEnv
``

Ensure that docker is installed
``
sudo apt get docker
``


Run the creation of the image named `rapgap` as:
``
sudo docker build -t rapgap .
``

Run the docker image as (`it` means interactively `rm` means that the container is removed after exiting:
``
sudo docker run -it --rm rapgap
``

Run rapgap inside the container
``
./testRun
``
