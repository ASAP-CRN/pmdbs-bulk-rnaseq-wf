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
		make \
		gcc \
		build-essential \
		libz-dev \
		libncurses5-dev \
		libbz2-dev \
		liblzma-dev \
	&& rm -rf /var/lib/apt/lists/*

# STAR
ARG STAR_VERSION
ENV STAR_VERSION "${STAR_VERSION}"
RUN wget https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz \
	&& tar -xzf ${STAR_VERSION}.tar.gz \
	&& mv STAR-${STAR_VERSION} /opt/star
RUN chmod a+x /opt/star \
	&& ln -s /opt/star /usr/local/bin/

# samtools
ARG SAMTOOLS_VERSION
ENV SAMTOOLS_VERSION "${SAMTOOLS_VERSION}"
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
	tar -jxvf samtools-${SAMTOOLS_VERSION}.tar.bz2 --directory /opt && \
	rm samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN cd /opt/samtools-${SAMTOOLS_VERSION} && \
	./configure --prefix=/usr/local/ && \
	make && \
	make install

COPY --from=scripts /opt/scripts /opt/scripts
ENV PATH "${PATH}:/opt/scripts"
