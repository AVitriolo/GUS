import smtplib
import email.mime.text
import email.mime.multipart
import email.mime.base
from email import encoders
import os

def send_email(subject, body, attachments, recipients, sender, app_password):
    for recipient in recipients:
        msg = email.mime.multipart.MIMEMultipart()
        msg["Subject"] = subject
        msg["From"] = sender
        msg["To"] = recipient

        msg.attach(email.mime.text.MIMEText(body, "plain"))

        for filepath in attachments:
            if not os.path.isfile(filepath):
                continue

            with open(filepath, "rb") as f:
                part = email.mime.base.MIMEBase("application", "octet-stream")
                part.set_payload(f.read())

            encoders.encode_base64(part)
            part.add_header(
                "Content-Disposition",
                f'attachment; filename="{os.path.basename(filepath)}"'
            )
            msg.attach(part)

        with smtplib.SMTP("smtp.gmail.com", 587) as server:
            server.starttls()
            server.login(sender, app_password)
            server.send_message(msg)
