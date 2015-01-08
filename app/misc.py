import unicodedata as ud

def simplify_string(odd_string):
    keepcharacters = (' ', '.', '_')
    safe_string = "".join(c for c in odd_string if c.isalnum()
                          or c in keepcharacters)\
        .rstrip().replace(' ', '_')
    if safe_string:
        safe_string = safe_string.encode('latin-1', 'ignore')\
            .decode('latin-1')
        safe_string = ud.normalize("NFC", safe_string)
    return safe_string
