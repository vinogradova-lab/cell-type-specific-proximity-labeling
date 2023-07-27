FROM continuumio/miniconda3

WORKDIR /work_dir

COPY environment.yml .
RUN conda env create -f myenv.yml

SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]

COPY scripts .

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "myenv", "python", "main_script.py"]
