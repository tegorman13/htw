---
title: simple ai chat
---


```{python}
#!pip install -q simpleaichat

from simpleaichat import AIChat
from getpass import getpass

```

```{python}
api_key = getpass("OpenAI Key: ")
```



```{python}

model = "gpt-3.5-turbo"  # in production, may want to use model="gpt-4" if have access

system_optimized = """You are a R code example.
Follow ALL the following rules:
- ONLY EVER RESPOND WITH CODE IN R MARKDOWN BLOCKS, AND NOTHING ELSE
- NEVER include code comments.
- Always assume library(tidyverse) is already executed.
"""


new_params = {
    "temperature": 0.1,
    "stop": ["```{r} ", "```\n"]
}
ai_2 = AIChat(api_key="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", system=system_optimized, model=model, params=new_params)


```


```{python}
%%time
response = ai_2("is_palindrome")
print(response)
```

print(ai_2)
{
  "id": "28b80cae-b99e-4d3f-aa95-c5a3847815cd",
  "created_at": "2023-06-11T21:38:38.572158+00:00",
  "auth": {
    "api_key": "**********"
  },
  "model": "gpt-3.5-turbo",
  "system": "You are a R code example.\nFollow ALL the following rules:\n- ONLY EVER RESPOND WITH CODE IN R MARKDOWN BLOCKS, AND NOTHING ELSE\n- NEVER include code comments.\n- Always assume library(tidyverse) is already executed.\n",
  "params": {
    "temperature": 0.1,
    "stop": [
      "```{r} ",
      "```\n"
    ]
  },
  "messages": [],
  "input_fields": [
    "content",
    "role"
  ],
  "save_messages": true,
  "total_prompt_length": 0,
  "total_completion_length": 0,
  "total_length": 0
}


```{python}
#print value of id in ai_2


```

```{python}



```