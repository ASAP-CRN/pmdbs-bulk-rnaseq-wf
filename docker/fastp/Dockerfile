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
		unzip \
		time \
	&& rm -rf /var/lib/apt/lists/*

# FastP
ARG FASTP_VERSION
ENV FASTP_VERSION "${FASTP_VERSION}"
RUN wget http://opengene.org/fastp/fastp.${FASTP_VERSION} \
	&& mv fastp.${FASTP_VERSION} /opt/fastp
RUN chmod a+x /opt/fastp \
	&& ln -s /opt/fastp /usr/local/bin/
