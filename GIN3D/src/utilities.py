import gi
import os
import traceback

from gi.repository import Gtk
gi.require_version('Gtk', '3.0')


def log(msg, mtype='e'):
    message_type = 'INFO' if mtype == 'i' else 'DEBUG' if mtype == 'd' else 'ERROR'
    print '[' + message_type + '] - ', msg
    if mtype == 'e':
        traceback.print_exc()


def display(msg, mtype, win):
    if mtype == 'e':
        dialog = Gtk.MessageDialog(win, 0, Gtk.MessageType.ERROR,
                                   Gtk.ButtonsType.CLOSE, "Oops! Error ...")
        dialog.format_secondary_text(msg)
        dialog.run()
        dialog.destroy()
    elif mtype == 'r':
        dialog = Gtk.MessageDialog(win, 0, Gtk.MessageType.ERROR,
                                   Gtk.ButtonsType.CLOSE, "Required value(s) missing ...")
        dialog.format_secondary_text(msg)
        dialog.run()
        dialog.destroy()
	  

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)
  
#def resource_path(relative):
#    return os.path.join(os.environ.get("_MEIPASS2", os.path.abspath(".")),relative)
	
def file_saveas(win, filters=None):
    dialog = Gtk.FileChooserDialog('Save As...', win, Gtk.FileChooserAction.SAVE,
                                   (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                                    Gtk.STOCK_SAVE_AS, Gtk.ResponseType.OK))
    dialog.set_do_overwrite_confirmation(True)

    if filters is not None:
        add_filters(dialog, filters)

    response = dialog.run()

    filename = None
    if response == Gtk.ResponseType.OK:
        filename = dialog.get_filename()
        log('File selected: ' + filename, 'd')
    elif response == Gtk.ResponseType.CANCEL:
        log('File selection canceled...', 'd')

    dialog.destroy()

    return filename


def file_chooser(win, filters=None):
    dialog = Gtk.FileChooserDialog('Open...', win, Gtk.FileChooserAction.OPEN,
                                   (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                                    Gtk.STOCK_OPEN, Gtk.ResponseType.OK))

    if filters is not None:
        add_filters(dialog, filters)

    response = dialog.run()
    filename = None
    if response == Gtk.ResponseType.OK:
        filename = dialog.get_filename()
        log('File selected: ' + filename, 'd')
    elif response == Gtk.ResponseType.CANCEL:
        log('File selection canceled...', 'd')

    dialog.destroy()

    return filename


def folder_chooser(win):
    dialog = Gtk.FileChooserDialog('Select...', win, Gtk.FileChooserAction.SELECT_FOLDER,
                                   (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                                    "Select", Gtk.ResponseType.OK))

    response = dialog.run()
    path = None
    if response == Gtk.ResponseType.OK:
        path = dialog.get_filename()
        log('Folder Selected: ' + str(path), 'd')
    elif response == Gtk.ResponseType.CANCEL:
        log('Folder selection canceled...', 'd')

    dialog.destroy()

    return path


def add_filters(dialog, filters):
    for k in filters:
        f = Gtk.FileFilter()
        f.set_name(k)
        f.add_pattern(filters[k])
        dialog.add_filter(f)


def first_child(p_id, c_id, index, gtkBuilder):
    p = gtkBuilder.get_object(p_id)
    grid_list = p.get_children()
    for c in grid_list[index + 1].get_children():
        if Gtk.Buildable.get_name(c) is not None:
            if Gtk.Buildable.get_name(c) == c_id + str(index + 1):
                return c
