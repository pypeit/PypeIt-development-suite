---
name: pyqt-expert
description: Master PyQt6 and PySide6 development for high-performance desktop applications. Handles native UI aesthetics, custom widgets, model-view architecture, multithreading with QThread, and application distribution. Use PROACTIVELY for desktop GUI design, system integration, or cross-platform desktop feature implementation.
---

## Use this skill when

- Developing desktop applications using Python (PyQt6 or PySide6).
- Designing native-looking UIs with custom styling (QSS) or Material Design.
- Implementing complex GUI patterns like Model/View (QAbstractItemModel).
- Handling long-running tasks without blocking the UI using QThread/QRunnable.
- Packaging and distributing Python desktop apps (PyInstaller, fbs, Nuitka).

## Instructions

- Prefer **PyQt6** or **PySide6** over older versions unless specified.
- Use **Model/View** architecture for large datasets or complex lists.
- Always separate UI logic from business logic.
- Use **QSS (Qt Style Sheets)** for modern, premium aesthetics.
- Implement **responsive layouts** using QLayout (QHBoxLayout, QVBoxLayout, QGridLayout).

## Capabilities

### UI Design & Aesthetics
- **Native Look & Feel**: Adapting to Windows, macOS, and Linux conventions.
- **Custom Styling**: Advanced QSS implementation, glassmorphism, and dark mode support.
- **Modern Themes**: Integration with `qt-material` or custom SVG-based icons.
- **Animations**: Using `QPropertyAnimation` and `QParallelAnimationGroup` for smooth transitions.

### Core Architecture
- **Signals & Slots**: Thread-safe communication between components.
- **Model/View Framework**: Custom models for `QTableView`, `QTreeView`, and `QListView`.
- **Multithreading**: Proper use of `QThread`, `QThreadPool`, and `signals` to keep UI responsive.
- **Event Loop Mastery**: Custom event filters and overriding `event()` for low-level control.

### System Integration
- **System Tray**: Implementing background apps with `QSystemTrayIcon`.
- **Native Dialogs**: Using `QFileDialog`, `QColorDialog`, and `QMessageBox` for system consistency.
- **Settings Persistence**: Managing app state with `QSettings` (INI, Registry, Plist).
- **Clipboard & Drag-and-Drop**: Managing deep system interactions.

### Deployment & Hardening
- **Packaging**: Freezing apps with `PyInstaller`, `briefcase`, or `Nuitka`.
- **Resource Management**: Compiling icons and images into `.qrc` / `_rc.py` files.
- **Hardening**: Integrating with `pyarmor` for code protection.

## Example Interactions
- "Create a PyQt6 dashboard with a side panel and glassmorphic tables."
- "Implement a background file processor in PyQt using QThread and a progress bar."
- "Style my PyQt application with a modern dark theme and custom SVG icons."
- "Set up a Model/View architecture for a searchable contact list."
