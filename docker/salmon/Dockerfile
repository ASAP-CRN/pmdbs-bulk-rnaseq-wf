FROM us-central1-docker.pkg.dev/dnastack-asap-parkinsons/workflow-images/util:1.1.1 as scripts

MAINTAINER Karen Fang <karen@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG GCLOUD_CLI_VERSION
ENV GCLOUD_CLI_VERSION "${GCLOUD_CLI_VERSION}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		wget \
		time \
		gzip \
		unzip \
		cmake \
		make \
		gcc \
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
	&& rm -rf /var/lib/apt/lists/*

# SALMON
ARG SALMON_VERSION
ENV SALMON_VERSION "${SALMON_VERSION}"
RUN wget https://github.com/COMBINE-lab/salmon/archive/v${SALMON_VERSION}.tar.gz \
	&& tar -xzf v${SALMON_VERSION}.tar.gz \
	&& mv salmon-${SALMON_VERSION} /opt/salmon
RUN cd /opt/salmon \
	&& cmake -DCMAKE_INSTALL_PREFIX=/usr/local/salmon \
	&& make \
	&& make install
