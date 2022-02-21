FROM python:3.9.9

WORKDIR .

RUN pip install --upgrade pip
COPY Project_4 .
RUN pip install -e project_spectra_temp/ flask
RUN pip install flask_sqlalchemy

LABEL "author"="Group_5" "email"="wdanqi@live.com"
LABEL build_date="2022-02-06"

#FLASK_PORT 5005
ENV FLASK_DEBUG=true FLASK_APP=run.py FLASK_ENV=development
ENV FLASK_PORT=5005
EXPOSE $FLASK_PORT

ENTRYPOINT [ "python", "frontend/run.py" ]




