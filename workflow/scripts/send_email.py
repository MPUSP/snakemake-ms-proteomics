#!/usr/bin/python3

# SEND EMAILS
# -----------------------------------------------------------------------------
#
# The purpose of this script is to send out reports
# that summarize the results of the pipeline by email.

from sys import exc_info
from os import path
from datetime import datetime
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders
from smtplib import SMTP


def text_body(snakemake):
    config_database = snakemake.params["config_database"]
    config_workflow = snakemake.params["config_workflow"]
    config_samplesheet = snakemake.params["config_samplesheet"]
    config_out_dir = path.abspath(snakemake.output["path"]).replace("/email", "")
    curr_time = datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S")
    report_html = snakemake.input["html"]
    report_pdf = snakemake.input["pdf"]
    protein = snakemake.input["protein"]
    comparison = snakemake.input["comparison"]
    text = f"""
        Hello,

        New proteomics data was processed with an automated analysis pipeline.
        The pipeline is based on fragpipe and MSstats, and can handle different
        types of MS data. The following input files were used:

          - sample sheet ('manifest'): {config_samplesheet}
          - workflow: {config_workflow}
          - database (organism's proteome): {config_database}

        The pipeline finished running at: {curr_time}

        The output was saved at this location: {config_out_dir}

        Detailed reports in PDF and HTML format are attached to this email.

        Bugs should be reported at: https://github.com/MPUSP/snakemake-ms-proteomics.
        The author(s) can be contacted via email at: jahn@mpusp.mpg.de.
        The pipeline is maintained by: Max-Planck-Unit for the Science of Pathogens, Berlin, Germany.

        Best regards!


        -------------------------------------------------------------------------
        Attachments:
        (Note: it is recommended to download the reports before opening)
        
        + Report in HTML format: {report_html}
        + Report in PDF format: {report_pdf}
        + Protein quantification table: {protein}
        + Comparison result table: {comparison}
        -------------------------------------------------------------------------
        """
    return text


def send_email(snakemake):
    # construct email
    email_body = text_body(snakemake)
    config = snakemake.params["config_email"]
    output_log = snakemake.log["path"]
    email_port = config.get("port")
    from_sender = config.get("from")
    msg = MIMEMultipart()
    msg["From"] = from_sender
    msg["To"] = ",".join(config.get("to"))
    msg["Subject"] = config.get("subject")
    msg.attach(MIMEText(email_body, "plain"))
    msg_noattach = msg.as_string()
    log = []

    # add reports as attachment to email
    for report in ["html", "pdf", "protein", "comparison"]:
        part = MIMEBase("application", "octet-stream")
        with open(snakemake.input[report], "rb") as file:
            part.set_payload(file.read())
        encoders.encode_base64(part)
        part.add_header(
            "Content-Disposition",
            "attachment; filename={}".format(path.abspath(snakemake.input[report])),
        )
        msg.attach(part)

    # send email
    if config.get("send"):
        try:
            server = SMTP(config.get("smtp_server"), email_port)
            server.starttls()
            server.login(config.get("smtp_user"), config.get("smtp_pw"))
            server.sendmail(from_sender, config.get("to"), msg.as_string())
            log += ["EMAIL: Report sent to: " + ", ".join(config.get("to"))]
            log += ["EMAIL: Module finished successfully"]
            error = 0
        except:
            error = str(exc_info())
            log += [
                "EMAIL, ERROR: " + error,
                "EMAIL: Report could not be sent. Check connection to email server.",
            ]
    else:
        error = 1
        log += [
            "EMAIL: Report has not been sent by mail. "
            + "Set 'email: send: True' in config.yml if you like to do so."
        ]

    # write log file
    with open(path.abspath(output_log), "w") as log_file:
        log_file.write("\n".join(log))
        if not error:
            log_file.write("\n\n +++ EMAIL BODY +++ \n\n")
            log_file.write(msg_noattach)

    return True


send_email(snakemake)
