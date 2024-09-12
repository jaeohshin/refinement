import os.path
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
from google.oauth2.credentials import Credentials
from googleapiclient import errors
from email.message import EmailMessage
import base64

__all__ = ['SubmitWithGmail']


class SubmitWithGmail:
    def __init__(self, **kwargs):
        """
        message: str
        subject: str
        to: str
        team_name: str
        token_path: path
        cred_path: path
        """
        self.__dict__.update(kwargs)

    def gmail_authenticate(self):
        SCOPES = ['https://mail.google.com/']
        creds = None
        if os.path.exists('token.json'):
            creds = Credentials.from_authorized_user_file(self.token_path, SCOPES)
        if not creds or not creds.valid:
            if creds and creds.expired and creds.refresh_token:
                creds.refresh(Request())
            else:
                flow = InstalledAppFlow.from_client_secrets_file(self.cred_path, SCOPES)
                creds = flow.run_local_server(port=0)
            with open('token.json', 'w') as token:
                token.write(creds.to_json())
        return build('gmail', 'v1', credentials=creds)

    @staticmethod
    def create_message(sender, to, subject, message_text):
        message = EmailMessage()
        message["From"] = sender
        message["To"] = to.split(",")
        message["Subject"] = subject
        message.set_content(message_text)
        return {'raw': base64.urlsafe_b64encode(message.as_bytes()).decode('utf8')}

    @staticmethod
    def send_message(service, user_id, message):
        try:
            message = service.users().messages().send(userId=user_id, body=message).execute()
            print('Message Id: %s' % message['id'])
            return message
        except errors.HttpError as error:
            print('An error occurred: %s' % error)

    def submit(self, title: str, message: str):
        service = self.gmail_authenticate()
        message = self.create_message(self.team_name, self.to, self.subject, self.message)
        self.send_message(service, "me", message)
