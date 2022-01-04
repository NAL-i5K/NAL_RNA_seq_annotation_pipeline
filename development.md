# The process of creating a new container for NAL RNA-Seq Annotation Pipeline on SCINet

A documentation on how to develop NAL RNA-Seq Annotation Pipeline. Since associating GitHub and Docker Hub accounts becomes a pro feature, the update process requires more steps as the following description.

## Prerequisite
- Docker

## Step 1: Build a Docker image on local

After merging the development branch into master on GitHub, build a Docker image from the GitHub repository on your local.
You can use the URL to build an image for a repository containing Dockerfile. [For detail](https://docs.docker.com/engine/reference/commandline/build/)
``` 
$  docker build github.com/NAL-i5K/NAL_RNA_seq_annotation_pipeline
``` 
Then rename the image and specify a tag.
``` 
$  docker image ls
$  docker tag [Image ID] [USERNAME/NAME-OF-REPO:TAG]  # the default tag is latest
``` 

## Step 2: Create a repository on Docker Hub and push the image to it

Log in to your Docker Hub account and create a repository (called NAME-OF-REPO). Then push the image built at step 1 to the repository on Docker Hub.
``` 
$  docker push [USERNAME/NAME-OF-REPO:TAG]
``` 
Check if the image is successfully uploaded to the Docker Hub repository.

## Step 3: Create a new container on SCINet

Connect to SCINet and navigate to the shared program directory. Import Docker images from Docker Hub to create a Singularity image and the corresponding container. [For detail](https://scinet.usda.gov/guide/singularity#4-docker-images)
``` 
$  singularity pull docker://[USERNAME/NAME-OF-REPO:TAG]
``` 

Congrats!!!