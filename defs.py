def groop(list_of,id):
    """
    Функция выбирает все записи с таким айди
    """
    gr = []
    for record in list_of:
        if record.id == id:
            gr.append(record)
    return gr
