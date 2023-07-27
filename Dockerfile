FROM continuumio/miniconda3

WORKDIR /work_dir

COPY environment_slim.yml .
RUN conda env create -f environment_slim.yml

SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]

COPY scripts .

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "myenv", "python", "main_script.py"]
