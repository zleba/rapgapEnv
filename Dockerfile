FROM ubuntu:20.04
RUN apt-get update
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install  wget vim g++ gfortran python make cmake texlive   ghostscript


RUN mkdir -p programs/lhapdfInstall programs/hepmcInstall  programs/rapgapInstall  programs/pythiaInstall programs/rapgapInstall

WORKDIR /programs

# Install LHAPDF

RUN wget -O  lhapdf-5.9.1.tar.gz https://lhapdf.hepforge.org/downloads?f=old/lhapdf-5.9.1.tar.gz
RUN tar -zxvf  lhapdf-5.9.1.tar.gz  && rm lhapdf-5.9.1.tar.gz

WORKDIR /programs/lhapdf-5.9.1

RUN  ./configure FCFLAGS="-std=legacy" --disable-pyext --prefix=`readlink -e ../lhapdfInstall` && make -j4 && make install


# Install HEPMC

WORKDIR /programs

RUN wget http://hepmc.web.cern.ch/hepmc/releases/hepmc2.06.11.tgz
RUN tar xvf hepmc2.06.11.tgz && rm hepmc2.06.11.tgz
WORKDIR /programs/HepMC-2.06.11
RUN ./configure --with-momentum=GEV --with-length=CM --prefix=/programs/hepmcInstall  && make -j4 && make install


# install Pythia

WORKDIR /programs

# The pythia6 originally from /afs/desy.de/group/alliance/mcg/public/tmp/pythia6/427.2
ADD pythia6 /programs/pythia6

WORKDIR /programs/pythia6
RUN ./configure   --prefix=/programs/pythiaInstall && make clean && make -j4 &&  make install 



# Install Rapgap

WORKDIR /programs


RUN wget --no-check-certificate  -O rapgap.tar.gz https://rapgap.hepforge.org/downloads/?f=rapgap-3.303.tar.gz
RUN tar -zxvf rapgap.tar.gz && rm rapgap.tar.gz

WORKDIR /programs/rapgap-3.303
RUN ./configure   --disable-shared --prefix=/programs/rapgapInstall --with-pythia6=/programs/pythiaInstall  --with-lhapdf=/programs/lhapdfInstall --with-hepmc=/programs/hepmcInstall

RUN make -j4
RUN make install


WORKDIR /programs/lhapdfInstall/share/lhapdf
RUN mkdir PDFsets && wget -O PDFsets/cteq6l.LHpdf https://lhapdf.hepforge.org/downloads?f=pdfsets/5.9.1//cteq6l.LHpdf

WORKDIR /programs/rapgapInstall/bin
ADD testRun /programs/rapgapInstall/bin
