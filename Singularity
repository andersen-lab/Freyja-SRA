Bootstrap: docker
From: mambaorg/micromamba:0.25.1

%labels
    Author Dylan Pilz
    Email dpilz@scripps.edu

%environment
    PATH /opt/conda/envs/freyja-sra/bin:$PATH

%post
    # Create the environment
    micromamba create -n freyja-sra

    # Install packages from environment.yml
    #cp $MAMBA_USER/environment.yml /tmp/environment.yml
    micromamba install -y -n freyja-sra -f /tmp/environment.yml && \
    micromamba clean --all --yes

%files
    environment.yml /tmp/environment.yml
