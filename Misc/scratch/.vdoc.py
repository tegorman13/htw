# type: ignore
# flake8: noqa
#
#
#
#
#
#
#!pip install -q simpleaichat

from simpleaichat import AIChat
from getpass import getpass

#
#
#
api_key = getpass("OpenAI Key: ")
#
#
#
#

params = {"temperature": 0.0}  # for reproducibility
model = "gpt-3.5-turbo"  # in production, may want to use model="gpt-4" if have access

ai = AIChat(api_key=api_key, console=False, params=params, model=model)
response = ai("Write an is_palindrome() function in Python.")
print(response)

#
#
#

model = "gpt-3.5-turbo"  # in production, may want to use model="gpt-4" if have access

system_optimized = """You are a R code example.
Follow ALL the following rules:
- ONLY EVER RESPOND WITH CODE IN R MARKDOWN BLOCKS, AND NOTHING ELSE
- NEVER include code comments.
- Always assume library(tidyverse) is already executed.
"""


new_params = {
    "temperature": 0.0,
    "stop": ["```{r} ", "```\n"]
}
ai_2 = AIChat(api_key="response = ai("Write an is_palindrome() function in Python.")
print(response)", system=system_optimized, model=model, params=new_params)


#
#
#
#
%%time
response = ai_2("is_palindrome")
print(response)
#
#
#
