
FROM gitpod/workspace-full

RUN brew install R
RUN pip install rpy2 --no-binary :all:
