#-*- coding:utf-8 -*-

fontColors = {
    'black': 30,
    'red': 31,
    'green': 32,
    'yellow': 33,
    'blue': 34,
    'magenta': 35,
    'cyan': 36,
    'lightgray': 37,
    'default': 39,
    'darkgray': 90,
    'lightred': 91,
    'lightgreen': 92,
    'lightyellow': 93,
    'lightblue': 94,
    'lightmagenta': 95,
    'lightcyan': 96,
    'white': 97
}

backgroundColors = {
    'black': 40,
    'red': 41,
    'green': 42,
    'yellow': 43,
    'blue': 44,
    'magenta': 45,
    'cyan': 46,
    'lightgray': 47,
    'default': 49,
    'darkgray': 100,
    'lightred': 101,
    'lightgreen': 102,
    'lightyellow': 103,
    'lightblue': 104,
    'lightmagenta': 105,
    'lightcyan': 106,
    'white': 107
}

formatting = {
    'normal': 0,
    'bold': 1,
    'light': 2,
    'italic': 3,
    'underline': 4,
    'blink': 5,
    'reverse': 7
}

def format(text, fc='red', bc='default', ft='normal'):
    """
    fc: font color
    bc: background color
    ft: formatting
    
    Colors: 
        blue, lightgray, lightyellow, default, darkgray, yellow, lightred, lightcyan, black, lightmagenta, lightblue, cyan, green, magenta, lightgreen, white, red
    Formattings:
        bold, normal, light, blink, italic, underline, reverse
    """
    code = "%d;%d;%d" % (formatting[ft], fontColors[fc], backgroundColors[bc])
    return "\x1b["+code+"m"+text+"\x1b[0m"

def f(text, fc='red', bc='default', ft='normal'):
    """
    fc: font color
    bc: background color
    ft: formatting
    
    Colors: 
        blue, lightgray, lightyellow, default, darkgray, yellow, lightred, lightcyan, black, lightmagenta, lightblue, cyan, green, magenta, lightgreen, white, red
    Formattings:
        bold, normal, light, blink, italic, underline, reverse
    """
    return format(text, fc=fc, bc=bc, ft=ft)
