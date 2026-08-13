# TOOLS.md - Local Notes

Skills define _how_ tools work. This file is for _your_ specifics — the stuff that's unique to your setup: camera names and locations, SSH hosts and aliases, preferred TTS voices, speaker/room names, device nicknames, anything environment-specific.

## Examples

```markdown
### Cameras

- living-room → Main area, 180° wide angle
- front-door → Entrance, motion-triggered

### SSH

- home-server → 192.168.1.100, user: admin

### TTS

- Preferred voice: "Nova" (warm, slightly British)
- Default speaker: Kitchen HomePod
```

## Why Separate?

Skills are shared. Your setup is yours. Keeping them apart means you can update skills without losing your notes, and share skills without leaking your infrastructure.

---

Add whatever helps you do your job. This is your cheat sheet.

## Related

- [Agent workspace](/concepts/agent-workspace)

## Voice replies (Telegram)

Voice replies are automatic. When the incoming message is a voice note, OpenClaw
synthesises your final text reply to speech and sends it as a voice message on its own
(`messages.tts.auto = "inbound"`, Microsoft voice `en-GB-SoniaNeural`). Text messages get
text replies.

So: just answer normally. Do **not** call `sarvam_text_to_speech` and do not try to attach
an audio file in order to "reply with voice" — that duplicates the reply and fails. Never
tell the user you are unable to send audio; you already did.

Use `sarvam_text_to_speech` only when the user explicitly asks for generated audio — a
specific Indian language or Bulbul voice, or a file they want to keep. Write it under
`/home/azureuser/.openclaw/workspace/media/sarvam` (the tool does this by default); paths
outside the workspace media root are rejected on send.

### Sarvam tools

- `sarvam_speech_to_text` — transcription; also wired as the automatic transcription hook
- `sarvam_translate` / `sarvam_transliterate` / `sarvam_identify_language` — text tools
- Sarvam has no British or American voices; all Bulbul voices are Indian-English or
  Indian-language. The British voice comes from Microsoft, not Sarvam.
