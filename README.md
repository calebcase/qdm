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
./qdm
```

---

[qdm-docker-badge]: https://img.shields.io/docker/cloud/build/calebcase/qdm
[qdm-docker]: https://hub.docker.com/repository/docker/calebcase/qdm
