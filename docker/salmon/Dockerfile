# Dockerfile based on: https://github.com/COMBINE-lab/salmon/blob/v1.10.3/docker/Dockerfile

ARG UBUNTU_VERSION

FROM us-central1-docker.pkg.dev/dnastack-asap-parkinsons/workflow-images/util:1.2.0 as scripts

FROM ubuntu:${UBUNTU_VERSION} as base
ENV UBUNTU_VERSION "${UBUNTU_VERSION}"

LABEL MAINTAINER="Karen Fang <karen@dnastack.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

# Salmon binary will be installed in /home/salmon/bin/salmon
WORKDIR /home

RUN apt-get -qq update \
	&& apt remove -y libcurl4 \
	&& apt-get install -y --no-install-recommends \
		curl \
		wget \
		unzip \
		git \
		gcc \
		make \
		g++ \
		libboost-all-dev \
		liblzma-dev \
		libbz2-dev \
		zlib1g-dev \
		libcurl4-openssl-dev \
		autoconf \
		apt-transport-https \
		gnupg \
		software-properties-common \
	&& apt-get clean

# Set up
RUN wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | apt-key add -

RUN apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main' \
	&& apt-get -qq update \
	&& apt-key --keyring /etc/apt/trusted.gpg del C1F34CDD40CD72DA \
	&& apt-get install -y \
		kitware-archive-keyring \
		cmake

# SALMON
ARG SALMON_VERSION
ENV SALMON_VERSION "${SALMON_VERSION}"
RUN wget https://github.com/COMBINE-lab/salmon/archive/v${SALMON_VERSION}.tar.gz \
	&& tar -xzf v${SALMON_VERSION}.tar.gz
RUN cd salmon-${SALMON_VERSION} \
	&& mkdir build \
	&& cd build \
	&& cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local/salmon \
	&& make \
	&& make install

FROM ubuntu:${UBUNTU_VERSION}
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		libhwloc5 \
		xxd \
		time \
		wget \
		ca-certificates \
	&& rm -rf /var/lib/apt/lists/*

COPY --from=base /usr/local/salmon/ /usr/local/
ENV PATH "${PATH}:/home/salmon/bin"
ENV LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:/usr/local/lib"

RUN echo "export PATH=$PATH" > /etc/environment
RUN echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH" > /etc/environment

# gcloud sdk; needed to upload output files
ARG GCLOUD_CLI_VERSION
ENV GCLOUD_CLI_VERSION "${GCLOUD_CLI_VERSION}"
RUN wget "https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${GCLOUD_CLI_VERSION}-linux-x86_64.tar.gz" \
	&& tar -zxvf "google-cloud-cli-${GCLOUD_CLI_VERSION}-linux-x86_64.tar.gz" --directory /opt \
	&& rm "google-cloud-cli-${GCLOUD_CLI_VERSION}-linux-x86_64.tar.gz"

ENV PATH "${PATH}:/opt/google-cloud-sdk/bin"

ARG PYTHON3_VERSION
ENV PYTHON3_VERSION "${PYTHON3_VERSION}"
RUN curl -O https://www.python.org/ftp/python/${PYTHON3_VERSION}/Python-${PYTHON3_VERSION}.tar.xz \
	&& tar -xvf Python-${PYTHON3_VERSION}.tar.xz --directory /opt/ \
	&& rm Python-${PYTHON3_VERSION}.tar.xz
RUN cd /opt/Python-${PYTHON3_VERSION} \
	&& ./configure \
	&& make \
	&& make altinstall

RUN update-alternatives --install /usr/bin/python3 python3 /opt/Python-${PYTHON3_VERSION}/python 1 \
	&& update-alternatives --set python3 /opt/Python-${PYTHON3_VERSION}/python

COPY --from=scripts /opt/scripts /opt/scripts
ENV PATH "${PATH}:/opt/scripts"
