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

ARG GCLOUD_CLI_VERSION
ENV GCLOUD_CLI_VERSION "${GCLOUD_CLI_VERSION}"

# Salmon binary will be installed in /home/salmon/bin/salmon
WORKDIR /home

RUN apt-get -qq update \
	&& apt remove -y libcurl4 \
	&& apt-get install -y --no-install-recommends \
		time \
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
		ca-certificates \
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
	&& rm -rf /var/lib/apt/lists/*

COPY --from=base /usr/local/salmon/ /usr/local/
ENV PATH "${PATH}:/home/salmon/bin"
ENV LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:/usr/local/lib"

RUN echo "export PATH=$PATH" > /etc/environment
RUN echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH" > /etc/environment

COPY --from=scripts /usr/lib/google-cloud-sdk /opt/google-cloud-sdk
COPY --from=scripts /opt/scripts /opt/scripts
ENV PATH "${PATH}:/opt/google-cloud-sdk/bin:/opt/scripts"
