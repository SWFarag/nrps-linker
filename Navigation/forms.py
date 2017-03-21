from flask_wtf import Form
from wtforms import TextAreaField, SubmitField,StringField, validators
class ContactForm(Form):
  name = StringField("Name",  [validators.required("Please enter your name.")])
  email = StringField("Email",  [validators.required("Please enter your email address."), validators.Email("Please enter a valid email address")])
  subject = StringField("Subject",  [validators.required("Please enter a subject.")])
  message = TextAreaField("Message",  [validators.required("Please enter a message.")])
  submit = SubmitField("Send")