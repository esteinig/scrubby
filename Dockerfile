FROM nvidia/cuda:12.1.0-cudnn8-devel-ubuntu20.04

# Set environment variables for non-interactive installation
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# Install necessary dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    bzip2 \
    ca-certificates \
    git \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
    sudo \
    tzdata \
    build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Configure the time zone
RUN ln -fs /usr/share/zoneinfo/$TZ /etc/localtime && dpkg-reconfigure --frontend noninteractive tzdata

# Install Micromamba
RUN curl -L https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

# Set up Micromamba
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
RUN ./bin/micromamba shell init -s bash -p $MAMBA_ROOT_PREFIX
ENV PATH=$MAMBA_ROOT_PREFIX/bin:$PATH

# Create a new environment and install PyTorch and its dependencies
RUN micromamba create -n pytorch -c pytorch -c nvidia -c conda-forge pytorch=2.3.0 pytorch-cuda=12.1 gxx torchvision torchaudio rust cudatoolkit -y

# Create an entrypoint script to activate the environment
RUN echo '#!/bin/bash\nexport PATH=/opt/micromamba/envs/pytorch/bin:${PATH}\nexec bash -c "$*"' > /entrypoint.sh
RUN chmod +x /entrypoint.sh

# Set the entrypoint
ENTRYPOINT ["/entrypoint.sh"]

# Set the working directory
WORKDIR /workspace

# Copy your application code to the container
COPY . /workspace

COPY Cargo.toml Cargo.toml
COPY src/ src/

# Uses the Pytorch GPU 'libtorch' to compile 'tch-rs' build script with Pytorch GPU support:
ENV LIBTORCH_USE_PYTORCH=1
ENV LD_LIBRARY_PATH=/opt/micromamba/envs/pytorch/lib/python3.12/site-packages/torch/lib:\$LD_LIBRARY_PATH
ENV PATH=/opt/micromamba/envs/pytorch/bin:/workspace/target/release:$PATH

RUN cargo build --release --features nn

# COPY pyproject.toml pyproject.toml
# COPY scrubby/ scrubby/

# RUN pip install .

# Set the default command
CMD ["bash"]
