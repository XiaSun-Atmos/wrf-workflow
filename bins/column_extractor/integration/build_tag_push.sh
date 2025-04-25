###
### Build local version of column-extractor container for code version tag
### Tag that local container build with the code version tag
### Also tag with necessary secondary tags for pushing to Harbor + AWS
### Then push containers to each repository
###
### To build v1.1.0, tag it, and push it to harbor-dev and AWS ECR
### in GSL's dev2 environment:
###
### bash ./integration/build_tag_push.sh -v v1.1.0 -b yes -p gsl-dev2
###
### (include optional "-t <PRODUCTION,DEVELOPMENT>" as extra tags)
###

## get command line arguments

while getopts v:b:t:p: flag

do

    case "${flag}" in
        v) VERSION_TAG=${OPTARG};;   # which version tag to use (ie v1.0.3, etc.)
        b) BUILD=${OPTARG};;         # whether we need to build image locally first
        t) EXTRA_TAGS+=(${OPTARG});; # any extra tags to add (ie PRODUCTION, DEVELOPMENT, etc)
        p) AWS_PROFILE=${OPTARG};;   # which AWS profile to use (ie gsl-dev2)
    esac

done

CONTAINER_TAG="column_extractor:$VERSION_TAG"

echo "Version tag is:   $VERSION_TAG"
echo "Container tag is: $CONTAINER_TAG"
echo "--------------------------"
echo

## build if needed

if [ $BUILD ]; then

    echo "Local build needed..."

    echo "Checking out main branch..."

    # first, checkout main branch
    git checkout main

    echo "Creating branch build/$VERSION_TAG from $VERSION_TAG..."

    # check out requested tag to separate build branch
    git checkout tags/$VERSION_TAG -b build/$VERSION_TAG
    
    # build the image using this version of the code
    echo "Building $CONTAINER_TAG from $VERSION_TAG..."
    docker build --no-cache --pull --build-arg="ARG_CODE_VERSION=$VERSION_TAG" -t $CONTAINER_TAG .

    echo "Checking out main branch..."

    # switch back to main
    git checkout main

    echo "Deleting branch build/$VERSION_TAG..."

    # delete build branch
    git branch -D build/$VERSION_TAG

    echo

fi

# create array to hold all tags (version tag + any extra) for harbor-dev / AWS ECR
ALL_TAGS=("$VERSION_TAG")
# append extra tags
ALL_TAGS+=("${EXTRA_TAGS[@]}")

## tag and push to harbor-dev (GSL's container registry)

for TAG in "${ALL_TAGS[@]}"

do

    # construct harbor-dev tag
    HARBOR_TAG="harbor-dev.gsd.esrl.noaa.gov/dsg/column_extractor:$TAG"

    # tag for harbov-dev
    echo "Harbor-dev tag is: $HARBOR_TAG"
    docker tag $CONTAINER_TAG $HARBOR_TAG

    # push to harbor-dev
    echo "Pushing to Harbor-dev..."
    docker push $HARBOR_TAG

    echo

done

## tag and push to AWS Elastic Container Repository (where the container is accessible for ECS tasks)

# get AWS account ID / region for provided profile
ACCOUNT_ID=$(aws sts get-caller-identity --profile $AWS_PROFILE --query Account --output text)
REGION=$(aws configure get region --profile $AWS_PROFILE)

# authenticate Docker against ECR
echo "Authenticating Docker to AWS ECR..."
aws ecr get-login-password --profile $AWS_PROFILE | docker login --username AWS --password-stdin $ACCOUNT_ID.dkr.ecr.$REGION.amazonaws.com

for TAG in "${ALL_TAGS[@]}"

do

    # construct full AWS ECR tag
    AWS_ECR_TAG="$ACCOUNT_ID.dkr.ecr.$REGION.amazonaws.com/dsg-column-extractor:$TAG"

    # tag for AWS ECR
    echo "AWS ECR tag is: $AWS_ECR_TAG"
    docker tag $CONTAINER_TAG $AWS_ECR_TAG

    # push to AWS ECR
    echo "Pushing to AWS ECR..."
    docker push $AWS_ECR_TAG

done
