**Installation of Docker in Ubuntu**

# install (if needed)
sudo apt-get update
sudo apt-get install -y ca-certificates curl gnupg lsb-release
# Add Docker repo (official)
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] \
  https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io

# add your user to docker group (logout/login required)
sudo usermod -aG docker $USER

# test
docker run --rm hello-world

spiseq@ASUS-TUF-Gaming-F15-FX506HF-FX506HF:~$ groups
spiseq adm cdrom sudo dip plugdev users lpadmin docker
spiseq@ASUS-TUF-Gaming-F15-FX506HF-FX506HF:~$ docker run --rm hello-world
Unable to find image 'hello-world:latest' locally
latest: Pulling from library/hello-world
17eec7bbc9d7: Pull complete 
Digest: sha256:f7931603f70e13dbd844253370742c4fc4202d290c80442b2e68706d8f33ce26
Status: Downloaded newer image for hello-world:latest

Hello from Docker!
This message shows that your installation appears to be working correctly.

To generate this message, Docker took the following steps:
 1. The Docker client contacted the Docker daemon.
 2. The Docker daemon pulled the "hello-world" image from the Docker Hub.
    (amd64)
 3. The Docker daemon created a new container from that image which runs the
    executable that produces the output you are currently reading.
 4. The Docker daemon streamed that output to the Docker client, which sent it
    to your terminal.

To try something more ambitious, you can run an Ubuntu container with:
 $ docker run -it ubuntu bash

Share images, automate workflows, and more with a free Docker ID:
 https://hub.docker.com/

 ✔ Advantages (why most labs use Docker/Singularity now)

Perfect reproducibility
Once an image is built, the tool versions never change.

Zero conflicts
Each run uses an isolated filesystem + isolated conda env.

Your operating system is unaffected
Even if Ubuntu breaks packages, the container still works.

Shareable / portable
You can give your image to someone else → same pipeline works on their system.

Supports HPC with Singularity
Docker → Singularity conversion is standard for cluster use.

Easy to archive your workflow
Imagine rerunning the same analysis 3 years later — Docker is the only reliable method.

❌ Disadvantages

Building the image can take time.

Requires Docker knowledge.

Slight overhead in container startup (but negligible for real pipelines).

For more examples and ideas, visit:
 https://docs.docker.com/get-started/

spiseq@ASUS-TUF-Gaming-F15-FX506HF-FX506HF:~$ 
