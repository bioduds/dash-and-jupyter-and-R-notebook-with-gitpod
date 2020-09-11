
FROM gitpod/workspace-full

RUN brew install python R
RUN pip install rpy2 --no-binary :all:
RUN run.sh
