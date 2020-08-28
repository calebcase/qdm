[![Docker Status][qdm-docker-badge]][qdm-docker]

# QDM

## Building from Source

### Requirements

The following tools are required:

* git
* git-lfs
* docker

**NOTE** Git LFS is required. Please ensure it is installed *before* attempting
to clone the repository.

### Build

Build the binary using docker:

```bash
git clone https://github.com/calebcase/qdm.git
cd qdm
docker build -t qdm .
```

Optionally, extract the binary from the docker image:

```bash
docker create -it --name qdm qdm
docker cp qdm:/usr/local/bin/qdm qdm
docker rm -f qdm
```

## Running

Run the binary from within docker:

```bash
docker run -it --rm -v $(pwd):/data qdm
```

Run locally (if you have extracted the binary):

```bash
./qdm -i test/data/input.h5 -o test/data/output.h5
```

During the run the output file contains temporary data that is used for
checkpointing and the final phase of analysis. If the debug flag is set to
false, then this temporary data is removed at the end of processing.  However
this does not reclaim the space due to the way h5 files work. In order to
reclaim the unused space it is necessary to run h5repack:

```bash
# If you haven't installed the h5 tools yet:
apt-get install hdf5-tools

h5repack test/data/output.h5 test/data/output.h5.repack
```

This can reduce the output file size significantly.

---

[qdm-docker-badge]: https://img.shields.io/docker/cloud/build/calebcase/qdm
[qdm-docker]: https://hub.docker.com/repository/docker/calebcase/qdm
