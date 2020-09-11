
FROM gitpod/workspace-full

RUN brew install R
RUN conda install -c r rpy2