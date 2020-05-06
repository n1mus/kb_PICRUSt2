FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.



RUN conda create --yes --name picrust2 conda=4.8.3
RUN conda install --yes --name picrust2  --channel bioconda --channel conda-forge picrust2=2.3.0_b 

RUN pip install pandas numpy
RUN pip install dotmap
RUN pip install seaborn
RUN pip install fastcluster

# install plotly
RUN apt-get update
RUN apt-get install --yes gcc libgtk2.0-0 libgtk-3-0 libxss1 libasound2
RUN pip install --upgrade pip
RUN pip install plotly psutil requests
RUN conda install --yes --channel plotly plotly-orca
RUN apt-get install --yes xvfb

ENV PYTHONUNBUFFERED=True

RUN apt-get install --yes vim

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
