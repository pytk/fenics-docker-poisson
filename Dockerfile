FROM quay.io/fenicsproject/stable:2017.2.0
SHELL ["/bin/bash", "-c"]
RUN sudo pip3 install --upgrade pip && sudo pip3 install https://github.com/JuliaPy/pyjulia/archive/master.zip#egg=julia ptvsd
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.1/julia-1.1.0-linux-x86_64.tar.gz
RUN tar -xzf julia-1.1.0-linux-x86_64.tar.gz
RUN echo 'export PATH="/julia-1.1.0/bin:$PATH"' >> ~/.bashrc
ENV PATH /home/fenics/julia-1.1.0/bin:$PATH
RUN ["/bin/bash", "-c", "source ~/.bashrc"]
WORKDIR '/home/fenics/shared'
COPY initiate.jl .
CMD ["julia", "initiate.jl"]