
FROM gitpod/workspace-full

RUN brew install python R
RUN pip install rpy2 --no-binary :all:
RUN pip3 install dash
RUN pip3 install dash-daq
RUN pip3 install jupyter pandas matplotlib
RUN pip install dash-auth