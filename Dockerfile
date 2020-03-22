FROM ubuntu:19.10 AS builder

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -qy update

RUN apt-get -qy install \
      build-essential \
      cmake \
      curl \
      git \
      meson \
      ninja-build \
      pkg-config

RUN apt-get -qy install \
      libargtable2-dev \
      libhdf5-dev

# Build and install GSL

WORKDIR /src
RUN curl -L -o gsl-2.6.tar.gz http://ftpmirror.gnu.org/gsl/gsl-2.6.tar.gz \
 && tar -xf gsl-2.6.tar.gz

WORKDIR /src/gsl-2.6
RUN ./configure \
 && make \
 && make install \
 && ldconfig

# Build and install OSQP

WORKDIR /src/github.com/oxfordcontrol
RUN git clone \
      --depth=1 \
      --single-branch \
      --branch v0.6.0 \
      --recursive \
      https://github.com/oxfordcontrol/osqp

WORKDIR /src/github.com/oxfordcontrol/osqp
RUN mkdir build \
 && cd build \
 && cmake -G 'Unix Makefiles' .. -DDLONG=OFF \
 && cmake --build . \
 && cmake --build . --target install \
 && ldconfig
COPY osqp.pc /usr/local/lib/pkgconfig/osqp.pc

# Build and install inih

WORKDIR /src/github.com/benhoyt
RUN git clone \
      --depth=1 \
      --single-branch \
      --branch r48 \
      https://github.com/benhoyt/inih

WORKDIR /src/github.com/benhoyt/inih
RUN mkdir build \
 && meson build \
 && ninja -C build \
 && ninja -C build install \
 && ldconfig

# Build and install QDM

COPY . /src/github.com/calebcase/qdm
WORKDIR /src/github.com/calebcase/qdm

RUN rm -fr build \
 && mkdir build \
 && meson build \
      --prefix=/usr/local \
      --buildtype=release \
 && ninja -C build \
 && ninja -C build test \
 && ninja -C build install \
 && ldconfig

FROM scratch
WORKDIR /data
COPY --from=builder /usr/local/bin/qdm /usr/local/bin/qdm
ENTRYPOINT ["/usr/local/bin/qdm"]